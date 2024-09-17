#include "embtrefftz.hpp"
#include "monomialfespace.hpp"
#include <cfloat>
#include <meshaccess.hpp>

using namespace ngbla;
using namespace ngcomp;

/// @returns the number of threads, that nngsolve might use.
int getNumberOfThreads ()
{
  // if no task manager is available, ngsolve will run sequentially.
  return (task_manager) ? task_manager->GetNumThreads () : 1;
}

template <typename SCAL>
void reorderMatrixColumns (MatrixView<SCAL> &matrix,
                           const Array<DofId> &dof_nrs, LocalHeap &local_heap)
{
  const auto heap_reset = HeapReset (local_heap);

  const auto matrix_shape = matrix.Shape ();

  FlatArray<int> map (dof_nrs.Size (), local_heap);
  for (int i = 0; i < map.Size (); i++)
    map[i] = i;

  QuickSortI (dof_nrs, map);
  auto elmat_b2_copy = FlatMatrix<SCAL> (get<0> (matrix_shape),
                                         get<1> (matrix_shape), local_heap);
  elmat_b2_copy = matrix;
  for (auto i : Range (get<0> (matrix_shape)))
    matrix.Col (i) = elmat_b2_copy.Col (map[i]);
}

template <typename SCAL, typename TDIST>
inline void addIntegrationToElementMatrix (
    MatrixView<SCAL, RowMajor, size_t, size_t, TDIST> elmat,
    const Array<shared_ptr<BilinearFormIntegrator>> &bf_integrators,
    const MeshAccess &mesh_access, const ElementId &element_id,
    const FESpace &fes, const FESpace &test_fes, LocalHeap &local_heap)
{
  const HeapReset hr (local_heap);

  auto &trafo = mesh_access.GetTrafo (element_id, local_heap);

  auto &test_fel = test_fes.GetFE (element_id, local_heap);
  auto &trial_fel = fes.GetFE (element_id, local_heap);

  const bool mixed_mode = std::addressof (test_fes) != std::addressof (fes);

  bool symmetric_so_far = true;
  int bfi_ind = 0;
  while (bfi_ind < bf_integrators.Size ())
    {
      auto &bfi = bf_integrators[bfi_ind];
      bfi_ind++;
      if (bfi->DefinedOnElement (element_id.Nr ()))
        {
          auto &mapped_trafo = trafo.AddDeformation (
              bfi->GetDeformation ().get (), local_heap);
          try
            {
              if (mixed_mode)
                {
                  const auto &mixed_fel
                      = MixedFiniteElement (trial_fel, test_fel);
                  bfi->CalcElementMatrixAdd (mixed_fel, mapped_trafo, elmat,
                                             symmetric_so_far, local_heap);
                }
              else
                {
                  bfi->CalcElementMatrixAdd (test_fel, mapped_trafo, elmat,
                                             symmetric_so_far, local_heap);
                }
            }
          catch (ExceptionNOSIMD e)
            {
              elmat = 0.0;
              cout << IM (6) << "ExceptionNOSIMD " << e.What () << endl
                   << "switching to scalar evaluation" << endl;
              bfi->SetSimdEvaluate (false);
              bfi_ind = 0;
            }
        }
    }
}

template <typename SCAL>
void extractVisibleDofs (FlatMatrix<SCAL> &elmat, const ElementId &element_id,
                         const FESpace &fes, const FESpace &test_fes,
                         Array<DofId> &dofs, Array<DofId> &test_dofs,
                         LocalHeap &local_heap)
{
  Array<DofId> vdofs, vtest_dofs;
  fes.GetDofNrs (element_id, vdofs, VISIBLE_DOF);
  test_fes.GetDofNrs (element_id, vtest_dofs, VISIBLE_DOF);
  FlatMatrix<SCAL> velmat (vtest_dofs.Size (), vdofs.Size (), local_heap);
  for (int jj = 0; jj < dofs.Size (); jj++)
    for (int ii = 0; ii < test_dofs.Size (); ii++)
      {
        auto j = vdofs.Pos (dofs[jj]);
        auto i = vtest_dofs.Pos (test_dofs[ii]);
        if (i != size_t (-1) && j != size_t (-1))
          velmat (i, j) = elmat (ii, jj);
      }
  dofs = std::move (vdofs);
  test_dofs = std::move (vtest_dofs);
  elmat.Assign (velmat);
}

/// Fills the two creators with the sparsity pattern needed for
/// the Trefftz embedding.
template <typename SCAL, typename NZ_FUNC>
INLINE size_t fillTrefftzTableCreators (
    TableCreator<int> &creator, TableCreator<int> &creator2,
    const vector<optional<Matrix<SCAL>>> &ETmats, const MeshAccess &ma,
    const FESpace &fes, NZ_FUNC nz_from_elnr, const size_t offset)
{
  const size_t ndof = fes.GetNDof ();
  const size_t ne = ma.GetNE (VOL);
  // number of the next Trefftz dof to create
  size_t next_trefftz_dof = offset;
  for (auto ei : ma.Elements (VOL))
    {
      if (!ETmats[ei.Nr ()])
        continue;

      size_t nz = nz_from_elnr (ei.Nr ());
      Array<DofId> dnums;
      fes.GetDofNrs (ei, dnums, VISIBLE_DOF);
      bool hasregdof = false;
      for (DofId d : dnums)
        if (IsRegularDof (d))
          {
            creator.Add (ei.Nr (), d);
            hasregdof = true;
          }
      // assumption here: Either all or no dof is regular
      if (hasregdof)
        {
          for (int d = 0; d < nz; d++)
            creator2.Add (ei.Nr (), next_trefftz_dof++);
        }
    }

  for (int d = 0, hcnt = 0; d < ndof; d++)
    if (HIDDEN_DOF == fes.GetDofCouplingType (d))
      {
        creator.Add (ne + hcnt, d);
        creator2.Add (ne + hcnt++, next_trefftz_dof++);
      }
  return next_trefftz_dof - offset;
}

/// creates two tables, used in `Elmats2Sparse`,
/// with a block-diagonal pattern, i.e. no overlap of dofs between two
/// elements.
template <typename SCAL>
INLINE size_t
createTrefftzTables (Table<int> &table, Table<int> &table2,
                     const vector<optional<Matrix<SCAL>>> &ETmats,
                     shared_ptr<const FESpace> fes, const size_t hidden_dofs)
{
  const auto ma = fes->GetMeshAccess ();
  const size_t ne = ma->GetNE (VOL);
  const size_t dim = fes->GetDimension ();
  TableCreator<int> creator (ne + hidden_dofs);
  TableCreator<int> creator2 (ne + hidden_dofs);
  size_t prevdofs = 0;

  for (; !creator.Done (); creator++, creator2++)
    {
      prevdofs = fillTrefftzTableCreators (
          creator, creator2, ETmats, *ma, *fes,
          [&] (size_t el_nr) { return ETmats[el_nr]->Width (); }, 0);
    }
  table = creator.MoveTable ();
  table2 = creator2.MoveTable ();
  return prevdofs;
}

template <typename SCAL>
INLINE size_t createConstrainedTrefftzTables (
    Table<int> &table, Table<int> &table2,
    const vector<optional<Matrix<SCAL>>> &ETmats, const FESpace &fes,
    const FESpace &fes_constraint, const size_t ndof_trefftz,
    const size_t hidden_dofs)
{
  const auto ma = fes.GetMeshAccess ();
  const size_t ne = ma->GetNE (VOL);
  const size_t ndof_constraint = fes_constraint.GetNDof ();
  const size_t dim = fes.GetDimension ();
  const size_t constraint_ndof = fes_constraint.GetNDof ();
  TableCreator<int> creator (ne + hidden_dofs);
  TableCreator<int> creator2 (ne + hidden_dofs);
  size_t trefftz_ndof = 0;

  for (; !creator.Done (); creator++, creator2++)
    {
      // first compute the Trefftz dofs. The dof numbers of the Trefftz dofs
      // are shifted up by constraint_ndof, to avoid conflicts between Trefftz
      // and Constraint dofs.
      trefftz_ndof = fillTrefftzTableCreators (
          creator, creator2, ETmats, *ma, fes,
          [ndof_trefftz] (size_t _) { return ndof_trefftz; }, constraint_ndof);
      (*testout) << "created " << trefftz_ndof << " many trefftz dofs"
                 << std::endl;

      // then compute the Constraint dofs.
      for (auto ei : ma->Elements (VOL))
        {
          if (!ETmats[ei.Nr ()])
            continue;

          Array<DofId> constraint_dnums;
          fes_constraint.GetDofNrs (ei, constraint_dnums, VISIBLE_DOF);

          bool hasregdof = false;
          for (DofId d : constraint_dnums)
            if (IsRegularDof (d))
              {
                // creator.Add (ei.Nr (), d);
                hasregdof = true;
              }
          // assumption here: Either all or no dof is regular
          if (hasregdof)
            {
              for (DofId d : constraint_dnums)
                creator2.Add (ei.Nr (), d);
            }
        }

      for (int d = 0, hcnt = 0; d < ndof_constraint; d++)
        if (HIDDEN_DOF == fes_constraint.GetDofCouplingType (d))
          {
            creator.Add (ne + hcnt, d);
            creator2.Add (ne + hcnt++, d);
          }
    }

  table = creator.MoveTable ();
  table2 = creator2.MoveTable ();
  return trefftz_ndof + constraint_ndof;
}

template <typename SCAL>
INLINE void
fillSparseMatrixWithData (SparseMatrix<SCAL> &P,
                          const vector<optional<Matrix<SCAL>>> &ETmats,
                          const Table<int> &table, const Table<int> &table2,
                          const MeshAccess &ma, const size_t hidden_dofs)
{
  const size_t ne = ma.GetNE (VOL);
  P.SetZero ();
  for (auto ei : ma.Elements (VOL))
    if (ETmats[ei.Nr ()])
      {
        P.AddElementMatrix (table[ei.Nr ()], table2[ei.Nr ()],
                            *(ETmats[ei.Nr ()]));
      }

  SCAL one = 1;
  FlatMatrix<SCAL> I (1, 1, &one);
  for (int hd = 0; hd < hidden_dofs; hd++)
    P.AddElementMatrix (table[ne + hd], table2[ne + hd], I);
}

INLINE size_t countHiddenDofs (const FESpace &fes)
{
  size_t hidden_dofs = 0;
  const size_t ndof = fes.GetNDof ();
  for (DofId d : Range (ndof))
    // branchless version of
    // if (HIDDEN_DOF == fes.GetDofCouplingType (d))
    //   hidden_dofs++;
    hidden_dofs += (HIDDEN_DOF == fes.GetDofCouplingType (d));
  return hidden_dofs;
}

namespace ngcomp
{
  /// assembles a global sparse matrix from the given element matrices.
  /// @param ETMats vector of all element matrices
  /// @param fes non-Trefftz finite element space
  /// @tparam SCAL scalar type of the matrix entries
  template <typename SCAL>
  shared_ptr<BaseMatrix>
  Elmats2Sparse (const vector<optional<Matrix<SCAL>>> ETmats,
                 shared_ptr<const FESpace> fes)
  {
    const auto ma = fes->GetMeshAccess ();

    size_t hidden_dofs = countHiddenDofs (*fes);

    Table<int> table, table2;
    const size_t prevdofs
        = createTrefftzTables (table, table2, ETmats, fes, hidden_dofs);

    auto P = make_shared<SparseMatrix<SCAL>> (fes->GetNDof (), prevdofs, table,
                                              table2, false);
    fillSparseMatrixWithData (*P, ETmats, table, table2, *ma, hidden_dofs);

    return P;
  }

  /// assembles a global sparse matrix from the given element matrices.
  /// @param ETMats vector of all element matrices
  /// @param fes non-Trefftz finite element space
  /// @tparam SCAL scalar type of the matrix entries
  template <typename SCAL>
  shared_ptr<BaseMatrix>
  Elmats2SparseConstrainedTrefftz (const vector<optional<Matrix<SCAL>>> ETmats,
                                   shared_ptr<const FESpace> fes,
                                   shared_ptr<const FESpace> fes_constraint,
                                   size_t ndof_trefftz)
  {
    const auto ma = fes->GetMeshAccess ();

    size_t hidden_dofs = countHiddenDofs (*fes);

    Table<int> table, table2;
    const size_t constraint_plus_trefftz_dofs
        = createConstrainedTrefftzTables (table, table2, ETmats, *fes,
                                          *fes_constraint, ndof_trefftz,
                                          hidden_dofs);

    auto P = make_shared<SparseMatrix<SCAL>> (
        fes->GetNDof (), constraint_plus_trefftz_dofs, table, table2, false);
    fillSparseMatrixWithData (*P, ETmats, table, table2, *ma, hidden_dofs);

    return P;
  }
}

void calculateBilinearFormIntegrators (
    const SumOfIntegrals &bf,
    Array<shared_ptr<BilinearFormIntegrator>> bfis[4])
{
  for (auto icf : bf.icfs)
    {
      DifferentialSymbol &dx = icf->dx;
      bfis[dx.vb] += icf->MakeBilinearFormIntegrator ();
    }
}

void calculateLinearFormIntegrators (
    const SumOfIntegrals &lf, Array<shared_ptr<LinearFormIntegrator>> lfis[4])
{
  for (auto icf : lf.icfs)
    {
      DifferentialSymbol &dx = icf->dx;
      lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
    }
}

bool fesHasHiddenDofs (const FESpace &fes)
{
  const size_t ndof = fes.GetNDof ();
  for (DofId d = 0; d < ndof; d++)
    if (HIDDEN_DOF == fes.GetDofCouplingType (d))
      return true;
  return false;
}

bool bfIsDefinedOnElement (const SumOfIntegrals &bf,
                           const Ngs_Element &mesh_element)
{
  for (auto icf : bf.icfs)
    {
      if (icf->dx.vb == VOL)
        if ((!icf->dx.definedonelements)
            || (icf->dx.definedonelements->Test (mesh_element.Nr ())))
          return true;
    }
  return false;
}

namespace ngbla
{
#ifdef NGSTREFFTZ_USE_LAPACK
  void
  LapackSVD (SliceMatrix<double, ColMajor> A, SliceMatrix<double, ColMajor> U,
             SliceMatrix<double, ColMajor> V)
  {
    static Timer t ("LapackSVD");
    RegionTimer reg (t);
    ngbla::integer n = A.Width (), m = A.Height ();
    Vector<> S (min (n, m));
    Array<double> work (4 * m * m + 6 * m + m + 100);
    Array<int> iwork (max (n, m) * 9);
    ngbla::integer info;
    char jobz = 'A';
    ngbla::integer lda = A.Dist (), ldu = U.Dist (), ldv = V.Dist ();
    ngbla::integer lwork = work.Size ();

    // if(std::is_same<SCAL,double>::value)
    dgesdd_ (&jobz, &m, &n, A.Data (), &lda, S.Data (), U.Data (), &ldu,
             V.Data (), &ldv, work.Data (), &lwork, iwork.Data (), &info);
    if (info != 0)
      throw Exception ("something went wrong in the svd "
                       + std::to_string (info));
    A = 0.0;
    A.Diag (0) = S;
  }

  void LapackSVD (SliceMatrix<Complex, ColMajor> A,
                  SliceMatrix<Complex, ColMajor> U,
                  SliceMatrix<Complex, ColMajor> V)
  {
    static Timer t ("LapackSVD");
    RegionTimer reg (t);
    ngbla::integer n = A.Width (), m = A.Height ();
    Vector<> S (min (n, m));
    Array<Complex> work (4 * m * m + 6 * m + m + 100);
    Array<int> iwork (max (n, m) * 9);
    Array<double> rwork (5 * max (n, m) * min (n, m) + 5 * min (n, m));
    ngbla::integer info;
    char jobz = 'A';
    ngbla::integer lda = A.Dist (), ldu = U.Dist (), ldv = V.Dist ();
    ngbla::integer lwork = work.Size ();

    zgesdd_ (&jobz, &m, &n, A.Data (), &lda, S.Data (), U.Data (), &ldu,
             V.Data (), &ldv, work.Data (), &lwork, rwork.Data (),
             iwork.Data (), &info);
    if (info != 0)
      throw Exception ("something went wrong in the svd "
                       + std::to_string (info));
    A = 0.0;
    A.Diag (0) = S;
  }
#endif

  template <typename SCAL>
  void getSVD (FlatMatrix<SCAL> A, FlatMatrix<SCAL, ColMajor> U,
               FlatMatrix<SCAL, ColMajor> V)
  {
    Matrix<SCAL, ColMajor> AA = A;
    // Matrix<SCAL,ColMajor> AA(A.Height(),A.Width());
    // for(int i=0;i<A.Height();i++)
    // for(int j=0;j<A.Width();j++)
    // AA(i,j)= A(i,j);
#ifdef NGSTREFFTZ_USE_LAPACK
    LapackSVD (AA, U, V);
#else
    cout << "No Lapack, using CalcSVD" << endl;
    CalcSVD (AA, U, V);
#endif
    A = static_cast<SCAL> (0.0);
    // A.Diag(0)=AA.Diag();
    for (int i = 0; i < min (A.Width (), A.Height ()); i++)
      A (i, i) = AA (i, i);
  }
}

/// `A = U * Sigma * V`
/// @return pseudo inverse of A (as some `ngbla::Expr` type to avoid
/// allocations)
template <typename SCAL>
auto invertSVD (const FlatMatrix<SCAL, ColMajor> &U,
                const FlatMatrix<SCAL> &Sigma,
                const FlatMatrix<SCAL, ColMajor> &V, LocalHeap &lh,
                size_t nz = 0)
{
  // let Sigma have dimension (m, n)
  // then U has dimension (m, m),
  // and V has dimension (n, n)

  auto U_T = Trans (U);
  auto V_T = Trans (V);
  const auto [m, n] = Sigma.Shape ();
  FlatMatrix<SCAL> Sigma_inv (n, m, lh);
  Sigma_inv = 0.;
  auto Sigma_el = Sigma.Diag ().begin ();
  auto Sigma_inv_el = Sigma_inv.Diag ().begin ();

  // Sigma_inv.Diag () = 1.0 / elmat_a.Diag ();
  size_t i = 0;
  // while (Sigma_el != Sigma.Diag().end () && abs (*Sigma_el) > 1e-10)
  while (nz > 0 ? i < max (m, n) - nz : abs (*Sigma_el) > 1e-10)
    {
      *Sigma_inv_el = 1.0 / *Sigma_el;
      i++;
      Sigma_el++;
      Sigma_inv_el++;
    }
  return V_T * Sigma_inv * U_T;
}

/// calculates from the given space and linear form integrators a particular
/// solution. note: allocates the solution on the given local heap.
/// @param lfis arrays of linear form integrators, for `VOL`, `BND`, `BBND`,
/// `BBBND`
/// @return the element-local solution, allocated on the local heap `mlh`
/// @tparam T type of matrix expression
template <typename SCAL, typename T>
FlatVector<SCAL>
calculateParticularSolution (Array<shared_ptr<LinearFormIntegrator>> lfis[4],
                             const FESpace &test_fes, const ElementId ei,
                             const MeshAccess &ma, const Array<DofId> &dofs,
                             const size_t ndof_test,
                             const Expr<T> &inverse_elmat, LocalHeap &mlh)
{
  // shall not be deallocated after function scope ends
  FlatVector<SCAL> elsol (dofs.Size (), mlh);

  const HeapReset hr (mlh);
  const auto &test_fel = test_fes.GetFE (ei, mlh);
  const auto &trafo = ma.GetTrafo (ei, mlh);
  FlatVector<SCAL> elvec (ndof_test, mlh), elveci (ndof_test, mlh);
  elvec = static_cast<SCAL> (0.0);
  for (const auto vorb : { VOL, BND, BBND, BBBND })
    {
      for (const auto &lfi : lfis[vorb])
        {
          if (lfi->DefinedOnElement (ei.Nr ()))
            {
              auto &mapped_trafo
                  = trafo.AddDeformation (lfi->GetDeformation ().get (), mlh);
              lfi->CalcElementVector (test_fel, mapped_trafo, elveci, mlh);
              elvec += elveci;
            }
        }
    }
  elsol = inverse_elmat * elvec;
  return elsol;
}

namespace ngcomp
{
  mutex stats_mutex;

  template <typename SCAL>
  std::tuple<vector<optional<Matrix<SCAL>>>, shared_ptr<BaseVector>>
  EmbTrefftz (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
              shared_ptr<SumOfIntegrals> lf, double eps,
              shared_ptr<FESpace> test_fes, int tndof, bool getrange,
              std::map<std::string, Vector<SCAL>> *stats)
  {
    static Timer svdtt ("svdtrefftz");
    RegionTimer reg (svdtt);
    LocalHeap lh (getNumberOfThreads () * 1000 * 1000);

    if (eps == 0 && tndof == 0 && test_fes == nullptr)
      throw Exception ("Need to specify eps, tndof, or test_fes");

    bool mixed_mode = true;
    if (test_fes == nullptr)
      {
        mixed_mode = false;
        test_fes = fes;
      }

    auto ma = fes->GetMeshAccess ();
    size_t ne = ma->GetNE (VOL);
    size_t ndof = fes->GetNDof ();
    size_t dim = fes->GetDimension ();

    Array<shared_ptr<BilinearFormIntegrator>> bfis[4]; // VOL, BND, ...
    calculateBilinearFormIntegrators (*bf, bfis);

    Array<shared_ptr<LinearFormIntegrator>> lfis[4];
    if (lf)
      calculateLinearFormIntegrators (*lf, lfis);

    vector<optional<Matrix<SCAL>>> ETmats (ne);
    auto lfvec = make_shared<VVector<SCAL>> (ndof);
    lfvec->operator= (0.0);

    size_t active_elements = 0;
    size_t test_local_ndof = test_fes->GetNDof () / ne;
    Vector<SCAL> sing_val_avg;
    Vector<double> sing_val_max;
    Vector<double> sing_val_min;

    const bool has_hidden_dofs = fesHasHiddenDofs (*fes);

    ma->IterateElements (VOL, lh, [&] (Ngs_Element ei, LocalHeap &mlh) {
      // skip this element, if the bilinear form is not defined
      // on this element
      if (!bfIsDefinedOnElement (*bf, ei))
        return;

      Array<DofId> dofs, test_dofs;
      test_fes->GetDofNrs (ei, test_dofs);
      fes->GetDofNrs (ei, dofs);

      FlatMatrix<SCAL> elmat (test_dofs.Size (), dofs.Size (), mlh);

      elmat = 0.;
      addIntegrationToElementMatrix (elmat, bfis[VOL], *ma, ElementId (ei),
                                     *fes, *test_fes, mlh);

      if (has_hidden_dofs) // extract visible dofs, if needed
        extractVisibleDofs (elmat, ei, *fes, *test_fes, dofs, test_dofs, mlh);

      FlatMatrix<SCAL, ColMajor> U (test_dofs.Size (), mlh),
          Vt (dofs.Size (), mlh);
      getSVD<SCAL> (elmat, U, Vt);

      // either tndof is given, or try to determine local size of trefftz
      // space, must be at least dim(trefftz)>=dim(trial)-dim(test)
      int nz = 0;
      if (tndof)
        nz = tndof;
      else
        {
          nz = max<int> (dofs.Size () - test_dofs.Size (), 0);
          for (int i = 0; i < min (elmat.Width (), elmat.Height ()); i++)
            if (abs (elmat (i, i)) < eps)
              nz++;
        }
      if (dofs.Size () < nz)
        throw Exception ("tndof too large, nz: " + to_string (nz)
                         + ", dofs:" + to_string (dofs.Size ())
                         + ", test_dofs:" + to_string (test_dofs.Size ()));

      if (getrange)
        ETmats[ei.Nr ()]
            = make_optional<Matrix<SCAL>> (U.Cols (0, dofs.Size () - nz));
      else
        ETmats[ei.Nr ()] = make_optional<Matrix<SCAL>> (
            Trans (Vt.Rows (dofs.Size () - nz, dofs.Size ())));

      if (stats)
        {
          const lock_guard<mutex> lock (stats_mutex);
          if (sing_val_avg.Size () == 0)
            {
              sing_val_avg.SetSize (elmat.Height ());
              sing_val_max.SetSize (elmat.Height ());
              sing_val_min.SetSize (elmat.Height ());
              sing_val_avg = 0;
              sing_val_max = 0;
              sing_val_min = DBL_MAX;
            }
          active_elements += 1;
          for (size_t i = 0; i < elmat.Height (); i++)
            {
              sing_val_avg[i] += elmat (i, i);
              sing_val_max[i] = max (sing_val_max[i], abs (elmat (i, i)));
              sing_val_min[i] = min (sing_val_min[i], abs (elmat (i, i)));
            }
        }

      if (lf)
        {
          auto elsol = calculateParticularSolution<SCAL> (
              lfis, *test_fes, ei, *ma, dofs, test_dofs.Size (),
              invertSVD (U, elmat, Vt, mlh, nz), mlh);
          lfvec->SetIndirect (dofs, elsol);
        }
    });

    if (stats)
      {
        sing_val_avg /= active_elements;
        (*stats)["singavg"] = Vector<SCAL> (sing_val_avg);
        (*stats)["singmax"] = Vector<double> (sing_val_max);
        (*stats)["singmin"] = Vector<double> (sing_val_min);
      }

    return std::make_tuple (ETmats, lfvec);
  }

  template std::tuple<vector<optional<Matrix<double>>>, shared_ptr<BaseVector>>
  EmbTrefftz<double> (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
                      shared_ptr<SumOfIntegrals> lf, double eps,
                      shared_ptr<FESpace> test_fes, int tndof, bool getrange,
                      std::map<std::string, Vector<double>> *stats);
  template std::tuple<vector<optional<Matrix<Complex>>>,
                      shared_ptr<BaseVector>>
  EmbTrefftz<Complex> (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
                       shared_ptr<SumOfIntegrals> lf, double eps,
                       shared_ptr<FESpace> test_fes, int tndof, bool getrange,
                       std::map<std::string, Vector<Complex>> *stats);

  template <typename SCAL>
  tuple<vector<optional<Matrix<SCAL>>>, shared_ptr<ngla::BaseVector>>
  EmbTrefftz (const SumOfIntegrals &op, const FESpace &fes,
              const FESpace &fes_test, const ngfem::SumOfIntegrals &cop_lhs,
              const ngfem::SumOfIntegrals &cop_rhs,
              const FESpace &fes_constraint,
              shared_ptr<const ngfem::SumOfIntegrals> linear_form,
              const size_t ndof_trefftz)
  {
    auto mesh_access = fes.GetMeshAccess ();
    const size_t num_elements = mesh_access->GetNE (VOL);
    // #TODO what is a good size for the local heap?
    // For the moment: large enough constant size.
    LocalHeap local_heap = LocalHeap (getNumberOfThreads () * 1000 * 1000);
    const size_t dim = fes.GetDimension ();
    const size_t dim_constraint = fes_constraint.GetDimension ();

    // calculate the integrators for the three bilinear forms,
    // each for VOL, BND, BBND, BBBND, hence 4 arrays per bilnear form
    Array<shared_ptr<BilinearFormIntegrator>> op_integrators[4],
        cop_lhs_integrators[4], cop_rhs_integrators[4];
    calculateBilinearFormIntegrators (op, op_integrators);
    calculateBilinearFormIntegrators (cop_lhs, cop_lhs_integrators);
    calculateBilinearFormIntegrators (cop_rhs, cop_rhs_integrators);

    vector<optional<Matrix<SCAL>>> element_matrices (num_elements);

    const bool fes_has_hidden_dofs = fesHasHiddenDofs (fes);
    const bool fes_constraint_has_hidden_dofs
        = fesHasHiddenDofs (fes_constraint);

    (*testout) << "fes has ndof:" << fes.GetNDof () << std::endl;
    (*testout) << "fes_constraint has ndof:" << fes_constraint.GetNDof ()
               << std::endl;
    (*testout) << "trefftz space has ndof:" << ndof_trefftz * num_elements
               << std::endl;

    auto particular_solution_vec = make_shared<VVector<SCAL>> (fes.GetNDof ());
    particular_solution_vec->operator= (0.0);

    // solve the following linear system in an element-wise fashion:
    // L @ T1 = B for the unknown matrix T1,
    // with the given matrices:
    //     /   \    /   \
    //  A= |B_1| B= |B_2|
    //     | L |    | 0 |
    //     \   /    \   /
    mesh_access->IterateElements (
        VOL, local_heap,
        [&] (Ngs_Element mesh_element, LocalHeap &local_heap) {
          const ElementId element_id = ElementId (mesh_element);
          const FiniteElement &fes_element
              = fes.GetFE (element_id, local_heap);

          // skip this element, if the bilinear forms are not defined
          // on this element
          if (!bfIsDefinedOnElement (op, mesh_element)
              || !bfIsDefinedOnElement (cop_lhs, mesh_element)
              || !bfIsDefinedOnElement (cop_rhs, mesh_element))
            return;

          Array<DofId> dofs, dofs_test, dofs_constraint;
          fes.GetDofNrs (element_id, dofs);
          fes_test.GetDofNrs (element_id, dofs_test);
          fes_constraint.GetDofNrs (element_id, dofs_constraint);

          //     /   \    /   \
          //  A= |B_1| B= |B_2|
          //     | L |    | 0 |
          //     \   /    \   /
          // with B_1.shape == (ndof_constraint, ndof),
          // L.shape == (ndof, ndof)
          // thus A.shape == (ndof + ndof_constraint, ndof)
          const size_t ndof = dofs.Size ();
          const size_t ndof_test = dofs_test.Size ();
          const size_t ndof_constraint = dofs_constraint.Size ();
          auto elmat_a = FlatMatrix<SCAL> (ndof_test + ndof_constraint, ndof,
                                           local_heap);
          auto [elmat_b1, elmat_l] = elmat_a.SplitRows (ndof_constraint);

          //     /   \    /   \
          //  A= |B_1| B= |B_2|
          //     | L |    | 0 |
          //     \   /    \   /
          // with B_2.shape == (ndof_constraint, ndof_constraint),
          // and B.shape == ( ndof_constraint + ndof, ndof_constraint)
          auto elmat_b = FlatMatrix<SCAL> (ndof_test + ndof_constraint,
                                           ndof_constraint, local_heap

          );
          elmat_a = static_cast<SCAL> (0.);
          elmat_b = static_cast<SCAL> (0.);

          // elmat_b2 is a view into elamt_b
          MatrixView<SCAL> elmat_b2 = elmat_b.Rows (ndof_constraint);

          // the diff. operator L operates only on volume terms
          addIntegrationToElementMatrix (elmat_l, op_integrators[VOL],
                                         *mesh_access, element_id, fes,
                                         fes_test, local_heap);
          for (const auto vorb : { VOL, BND, BBND, BBBND })
            {
              addIntegrationToElementMatrix (
                  elmat_b1, cop_lhs_integrators[vorb], *mesh_access,
                  element_id, fes, fes_constraint, local_heap);
              addIntegrationToElementMatrix (
                  elmat_b2, cop_rhs_integrators[vorb], *mesh_access,
                  element_id, fes_constraint, fes_constraint, local_heap);
            }
          if (fes_has_hidden_dofs)
            throw std::invalid_argument (
                "fes has hidden dofs, not supported at the moment");
          if (fes_has_hidden_dofs)
            extractVisibleDofs (elmat_l, element_id, fes, fes, dofs, dofs,
                                local_heap);

          // reorder elmat_b2
          // #TODO is this really necessary?
          reorderMatrixColumns (elmat_b2, dofs_constraint, local_heap);

          //(*testout) << "elmat_a:\n" << elmat_a << std::endl;
          (*testout) << "calculated elmat_a on element:\n"
                     << element_id << std::endl;

          // singular value decomposition of elmat_a:
          // U * sigma * V = elmat_a
          // elmat_a gets overwritten with sigma
          // let elmat_a have dimension (m, n)
          // the U has dimension (m, m),
          // and V has dimension (n, n)
          FlatMatrix<SCAL, ColMajor> U (ndof_constraint + ndof_test,
                                        local_heap),
              V (ndof, local_heap);
          getSVD<SCAL> (elmat_a, U, V);

          // We use the inverse of `elmat_a` twice. In the calculation
          // of the particular solution, we only use the right block of the
          // matrix. Therefore, it is beneficial to store the whole matrix, not
          // just the multiplication expression.
          const auto elmat_a_inv_expr = invertSVD (U, elmat_a, V, local_heap);
          FlatMatrix<SCAL> elmat_a_inv (ndof, ndof_constraint + ndof_test,
                                        local_heap);
          // Calculate the matrix entries and write them to memory.
          elmat_a_inv = elmat_a_inv_expr;

          // P = (T1 | T2)
          Matrix<SCAL> elmat_p (ndof, ndof_trefftz + ndof_constraint);
          // T1 has dimension (ndof, ndof_constraint)
          // T2 has dimension (ndof, ndof)
          auto [elmat_t1, elmat_t2] = elmat_p.SplitCols (ndof_constraint);

          // T1 solves A @ T1 = B,
          // i.e. T1 = A^{-1} @ B.
          // A has dimension (ndof + ndof_constraint, ndof),
          // B has dimension (ndof + ndof_constraint, ndof_constraint),
          // so T1 has dimension (ndof, ndof_constraint)
          elmat_t1 = elmat_a_inv * elmat_b;
          //(*testout) << "t1 from element " << element_id << "\n"
          //           << elmat_t1 << std::endl;
          (*testout) << "calculated t1 from element " << element_id << "\n"
                     << std::endl;

          elmat_t2 = Trans (V.Rows (ndof_test - ndof_trefftz, ndof_test));
          //(*testout) << "t2 from element " << element_id << "\n"
          //           << elmat_t2 << std::endl;
          (*testout) << "calculated t2 from element " << element_id << "\n"
                     << std::endl;

          //(*testout) << "elmat_p shape: (" << ndof << ", "
          //           << ndof_constraint + ndof_trefftz << ")" << std::endl;

          //(*testout) << "p from element " << element_id << "\n"
          //           << elmat_p << std::endl;
          (*testout) << "calculated p from element " << element_id << "\n"
                     << std::endl;
          element_matrices[element_id.Nr ()]
              = make_optional<Matrix<SCAL>> (elmat_p);
          (*testout) << "allocated p as a shared pointer for element "
                     << element_id << "\n"
                     << std::endl;

          if (linear_form)
            {
              // the particular solution only needs to be computed for the
              // Trefftz part of the element diff. operator.
              const auto [_, elmat_l_inv]
                  = elmat_a_inv.SplitCols (ndof_constraint);
              (*testout) << "got elmat_l pseudoinverse for element "
                         << element_id << "\n"
                         << std::endl;

              Array<shared_ptr<LinearFormIntegrator>> lfis[4];
              calculateLinearFormIntegrators (*linear_form, lfis);
              const auto part_sol = calculateParticularSolution<SCAL> (
                  lfis, fes_test, element_id, *mesh_access, dofs, ndof_test,
                  elmat_l_inv, local_heap);
              (*testout) << "calculated particular solution on element "
                         << element_id << "\n"
                         << std::endl;
              particular_solution_vec->SetIndirect (dofs, part_sol);
              (*testout) << "wrote particular solution on element "
                         << element_id << "\n"
                         << std::endl;
            }
        });

    return make_tuple<> (element_matrices, particular_solution_vec);
  }

  ////////////////////////// EmbTrefftzFESpace ///////////////////////////

  template <typename T>
  shared_ptr<BaseVector>
  EmbTrefftzFESpace<T>::SetOp (shared_ptr<SumOfIntegrals> bf,
                               shared_ptr<SumOfIntegrals> lf, double eps,
                               shared_ptr<FESpace> test_fes, int tndof)
  {
    static Timer timer ("EmbTrefftz: SetOp");

    // shared_ptr<FESpace> use_as_l2 =
    // dynamic_pointer_cast<FESpace>(const_cast<EmbTrefftzFESpace*>(this)->shared_from_this());

    shared_ptr<BaseVector> lfvec;
    if (!this->IsComplex ())
      {
        auto embtr = EmbTrefftz<double> (bf, fes, lf, eps, test_fes, tndof,
                                         false, nullptr);
        this->ETmats = std::get<0> (embtr);
        lfvec = std::get<1> (embtr);
      }
    else
      {
        auto embtr = EmbTrefftz<Complex> (bf, fes, lf, eps, test_fes, tndof,
                                          false, nullptr);
        this->ETmatsC = std::get<0> (embtr);
        lfvec = std::get<1> (embtr);
      }

    adjustDofsAfterSetOp ();
    return lfvec;
  }

  template <typename T>
  shared_ptr<BaseVector>
  EmbTrefftzFESpace<T>::SetOp (shared_ptr<const SumOfIntegrals> op,
                               shared_ptr<const SumOfIntegrals> cop_lhs,
                               shared_ptr<const SumOfIntegrals> cop_rhs,
                               shared_ptr<const FESpace> fes_constraint,
                               shared_ptr<const FESpace> fes_test,
                               shared_ptr<const SumOfIntegrals> linear_form,
                               size_t ndof_trefftz)
  {
    static Timer timer ("EmbTrefftz: SetOp");

    shared_ptr<BaseVector> particular_solution = nullptr;

    if (!fes || !cop_lhs || !cop_rhs || !fes_constraint)
      throw std::invalid_argument ("All pointers except for op, fes_test and "
                                   "linear_form may not be null.");

    // fes_test may be null. If it is null, then choose the trial space as the
    // test space as well.
    const FESpace &fes_test_ref = (fes_test) ? *fes_test : *fes;

    // if op is null, set op to be an empty sum.
    const SumOfIntegrals op_default{};
    if (!op)
      ndof_trefftz = 0;
    const SumOfIntegrals &op_ref = (op) ? *op : op_default;

    if (!this->IsComplex ())
      {
        std::tie (this->ETmats, particular_solution) = EmbTrefftz<double> (
            op_ref, *fes, fes_test_ref, *cop_lhs, *cop_rhs, *fes_constraint,
            linear_form, ndof_trefftz);
      }
    else
      {
        std::tie (this->ETmatsC, particular_solution)
            = EmbTrefftz<Complex> (*op, *fes, fes_test_ref, *cop_lhs, *cop_rhs,
                                   *fes_constraint, linear_form, ndof_trefftz);
      }

    adjustDofsAfterSetOp ();

    return particular_solution;
  }

  template <typename T> void EmbTrefftzFESpace<T>::adjustDofsAfterSetOp ()
  {
    T::Update ();
    int ndof = fes->GetNDof ();
    all2comp.SetSize (ndof);
    all2comp = 0;

    for (auto ei : this->ma->Elements (VOL))
      {
        int nz = this->IsComplex () ? (ETmatsC[ei.Nr ()])->Width ()
                                    : (ETmats[ei.Nr ()])->Width ();
        Array<DofId> dofs;
        T::GetDofNrs (ei, dofs);
        for (int i = nz; i < dofs.Size (); i++)
          all2comp[dofs[i]] = NO_DOF_NR_CONDENSE;
        ////for(int i=0;i<dofs.Size()-nz;i++)
        // for(int i=nz;i<dofs.Size();i++)
        ////this->free_dofs->Clear(dofs[i]);
        ////this->SetDofCouplingType(dofs[i],UNUSED_DOF );
        // this->SetDofCouplingType(dofs[i],HIDDEN_DOF );
      }

    int newndof = 0;
    for (DofId &d : all2comp)
      if (d == 0)
        d = newndof++;

    this->SetNDof (newndof);

    // this->UpdateDofTables();
    // this->UpdateCouplingDofArray();

    this->ctofdof.SetSize (newndof);
    for (int i = 0; i < ndof; i++)
      if (all2comp[i] >= 0)
        this->ctofdof[all2comp[i]] = fes->GetDofCouplingType (i);
    // this->ctofdof[all2comp[i]] = T::GetDofCouplingType (i);

    T::FinalizeUpdate ();
  }

  template <typename T>
  void EmbTrefftzFESpace<T>::GetDofNrs (ElementId ei, Array<int> &dnums) const
  {
    T::GetDofNrs (ei, dnums);
    if (all2comp.Size () == fes->GetNDof ())
      for (DofId &d : dnums)
        if (IsRegularDof (d))
          d = all2comp[d];
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformMR (ElementId ei, SliceMatrix<double> mat,
                                      TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: MTransform");
    RegionTimer reg (timer);

    size_t nz = (ETmats[ei.Nr ()])->Width ();
    Matrix<double> temp_mat (mat.Height (), mat.Width ());

    if (type == TRANSFORM_MAT_LEFT)
      {
        temp_mat.Rows (0, nz) = Trans (*(ETmats[ei.Nr ()])) * mat;
        mat = temp_mat;
      }
    if (type == TRANSFORM_MAT_RIGHT)
      {
        temp_mat.Cols (0, nz) = mat * (*(ETmats[ei.Nr ()]));
        mat = temp_mat;
      }
    if (type == TRANSFORM_MAT_LEFT_RIGHT)
      {
        temp_mat.Cols (0, nz) = mat * (*(ETmats[ei.Nr ()]));
        mat.Cols (0, nz).Rows (0, nz) = Trans (*(ETmats[ei.Nr ()])) * temp_mat;
      }
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformMC (ElementId ei, SliceMatrix<Complex> mat,
                                      TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: MTransform");
    RegionTimer reg (timer);

    size_t nz = (ETmatsC[ei.Nr ()])->Width ();
    Matrix<Complex> temp_mat (mat.Height (), mat.Width ());

    if (type == TRANSFORM_MAT_LEFT)
      {
        temp_mat.Rows (0, nz) = Trans (*(ETmatsC[ei.Nr ()])) * mat;
        mat = temp_mat;
      }
    if (type == TRANSFORM_MAT_RIGHT)
      {
        temp_mat.Cols (0, nz) = mat * (*(ETmatsC[ei.Nr ()]));
        mat = temp_mat;
      }
    if (type == TRANSFORM_MAT_LEFT_RIGHT)
      {
        temp_mat.Cols (0, nz) = mat * (*(ETmatsC[ei.Nr ()]));
        mat.Cols (0, nz).Rows (0, nz)
            = Trans (*(ETmatsC[ei.Nr ()])) * temp_mat;
      }
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformVR (ElementId ei, SliceVector<double> vec,
                                      TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: VTransform");
    RegionTimer reg (timer);

    size_t nz = (ETmats[ei.Nr ()])->Width ();

    if (type == TRANSFORM_RHS)
      {
        Vector<double> new_vec (nz);
        new_vec = Trans (*(ETmats[ei.Nr ()])) * vec;
        vec = new_vec;
      }
    else if (type == TRANSFORM_SOL)
      {
        Vector<double> new_vec (vec.Size ());
        new_vec = (*(ETmats[ei.Nr ()])) * vec.Range (0, nz);
        vec = new_vec;
      }
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformVC (ElementId ei, SliceVector<Complex> vec,
                                      TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: VTransform");
    RegionTimer reg (timer);

    size_t nz = (ETmatsC[ei.Nr ()])->Width ();

    if (type == TRANSFORM_RHS)
      {
        Vector<Complex> new_vec (nz);
        new_vec = Trans (*(ETmatsC[ei.Nr ()])) * vec;
        vec = new_vec;
      }
    else if (type == TRANSFORM_SOL)
      {
        Vector<Complex> new_vec (vec.Size ());
        new_vec = (*(ETmatsC[ei.Nr ()])) * vec.Range (0, nz);
        vec = new_vec;
      }
  }

  template <typename T>
  shared_ptr<GridFunction>
  EmbTrefftzFESpace<T>::Embed (shared_ptr<GridFunction> tgfu)
  {
    LocalHeap lh (1000 * 1000 * 1000);
    Flags flags;

    auto tvec = tgfu->GetVectorPtr ();

    auto gfu = CreateGridFunction (this->fes, "pws", flags);
    gfu->Update ();
    auto vec = gfu->GetVectorPtr ();

    this->ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
      Array<DofId> dofs, tdofs;
      this->fes->GetDofNrs (ei, dofs);
      tdofs.SetSize (dofs.Size ());
      tdofs.SetSize0 ();

      for (DofId d : dofs)
        if (all2comp[d] >= 0)
          tdofs.Append (all2comp[d]);

      if (this->IsComplex ())
        {
          FlatVector<Complex> telvec (tdofs.Size (), mlh);
          tvec->GetIndirect (tdofs, telvec);
          FlatVector<Complex> elvec (dofs.Size (), mlh);
          elvec = (*(ETmatsC[ei.Nr ()])) * telvec;
          vec->SetIndirect (dofs, elvec);
        }
      else
        {
          FlatVector<> telvec (tdofs.Size (), mlh);
          tvec->GetIndirect (tdofs, telvec);
          FlatVector<> elvec (dofs.Size (), mlh);
          elvec = (*(ETmats[ei.Nr ()])) * telvec;
          vec->SetIndirect (dofs, elvec);
        }
    });
    return gfu;
  }

  template <typename T>
  shared_ptr<BaseMatrix> EmbTrefftzFESpace<T>::GetEmbedding ()
  {
    if (this->IsComplex ())
      return Elmats2Sparse<Complex> (ETmatsC, this->fes);
    else
      return Elmats2Sparse<double> (ETmats, this->fes);
  }

  // template class EmbTrefftzFESpace<L2HighOrderFESpace,
  // shared_ptr<L2HighOrderFESpace>>;
  // static RegisterFESpace<
  // EmbTrefftzFESpace<L2HighOrderFESpace, shared_ptr<L2HighOrderFESpace>>>
  // initembt ("L2EmbTrefftzFESpace");
  // template class EmbTrefftzFESpace<MonomialFESpace,
  // shared_ptr<MonomialFESpace>>;
  // static RegisterFESpace<
  // EmbTrefftzFESpace<MonomialFESpace, shared_ptr<MonomialFESpace>>>
  // initembt3 ("MonomialEmbTrefftzFESpace");
}

////////////////////////// python interface ///////////////////////////

#ifdef NGS_PYTHON
template <typename T> void ExportETSpace (py::module m, string label)
{
  auto pyspace
      = ngcomp::ExportFESpace<ngcomp::EmbTrefftzFESpace<T>> (m, label);

  pyspace.def (py::init ([pyspace] (shared_ptr<T> fes) {
                 py::list info;
                 auto ma = fes->GetMeshAccess ();
                 info.append (ma);
                 auto nfes = make_shared<ngcomp::EmbTrefftzFESpace<T>> (fes);
                 nfes->Update ();
                 nfes->FinalizeUpdate ();
                 connect_auto_update (nfes.get ());
                 return nfes;
               }),
               py::arg ("fes"));

  pyspace
      .def ("SetOp",
            static_cast<shared_ptr<BaseVector> (EmbTrefftzFESpace<T>::*) (
                shared_ptr<SumOfIntegrals>, shared_ptr<SumOfIntegrals>, double,
                shared_ptr<FESpace>, int)> (
                &ngcomp::EmbTrefftzFESpace<T>::SetOp),
            py::arg ("bf"), py::arg ("lf") = nullptr, py::arg ("eps") = 0,
            py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0)
      .def ("SetOp",
            static_cast<shared_ptr<BaseVector> (EmbTrefftzFESpace<T>::*) (
                shared_ptr<const SumOfIntegrals>,
                shared_ptr<const SumOfIntegrals>,
                shared_ptr<const SumOfIntegrals>, shared_ptr<const FESpace>,
                shared_ptr<const FESpace>, shared_ptr<const SumOfIntegrals>,
                const size_t)> (&ngcomp::EmbTrefftzFESpace<T>::SetOp),
            py::arg ("op"), py::arg ("cop_lhs"), py::arg ("cop_rhs"),
            py::arg ("fes_constraint"), py::arg ("fes_test") = nullptr,
            py::arg ("linear_form") = nullptr, py::arg ("ndof_trefftz") = 0)
      .def ("Embed", &ngcomp::EmbTrefftzFESpace<T>::Embed)
      .def ("GetEmbedding", &ngcomp::EmbTrefftzFESpace<T>::GetEmbedding);
}

/// call `EmbTrefftz` for the ConstrainedTrefftz procedure and pack the
/// resulting element matrices in a sparse matrix.
/// Assumption: all shared pointers come from python, so they *should* be safe
/// to dereference. Exception to that: the `linear_form` can be the `nullptr`.
tuple<shared_ptr<BaseMatrix>, shared_ptr<ngla::BaseVector>>
pythonConstrTrefftzWithLf (optional<const SumOfIntegrals> op_maybe,
                           shared_ptr<const FESpace> fes,
                           shared_ptr<const SumOfIntegrals> cop_lhs,
                           shared_ptr<const SumOfIntegrals> cop_rhs,
                           shared_ptr<const FESpace> fes_constraint,
                           shared_ptr<const SumOfIntegrals> linear_form,
                           size_t ndof_trefftz,
                           shared_ptr<const FESpace> fes_test_ptr)
{
  // guard against unwanted segfaults by checking that the pointers are not
  // null.
  if (!fes || !cop_lhs || !cop_rhs || !fes_constraint)
    throw std::invalid_argument (
        "Some arguments passed were None, which are required to be not-None.");

  // if fes_test is null, i.e. no test space is given, use fes as trial and
  // test space
  const FESpace &fes_test = (fes_test_ptr) ? *fes_test_ptr : *fes;

  // if op_maybe is empty, we need a default op to pass to the C++ code
  const SumOfIntegrals op_default{};
  if (!op_maybe)
    ndof_trefftz = 0;
  const SumOfIntegrals &op = (op_maybe) ? *op_maybe : op_default;

  auto [P, u_lf]
      = EmbTrefftz<double> (op, *fes, fes_test, *cop_lhs, *cop_rhs,
                            *fes_constraint, linear_form, ndof_trefftz);
  return std::make_tuple (ngcomp::Elmats2SparseConstrainedTrefftz<double> (
                              P, fes, fes_constraint, ndof_trefftz),
                          u_lf);
}

/// call `EmbTrefftz` for the ConstrainedTrefftz procedure and pack the
/// resulting element matrices in a sparse matrix.
shared_ptr<BaseMatrix>
pythonConstrTrefftz (optional<const SumOfIntegrals> op,
                     shared_ptr<const FESpace> fes,
                     shared_ptr<const SumOfIntegrals> cop_lhs,
                     shared_ptr<const SumOfIntegrals> cop_rhs,
                     shared_ptr<const FESpace> fes_constraint,
                     const size_t ndof_trefftz,
                     shared_ptr<const FESpace> fes_test)
{
  return std::get<0> (pythonConstrTrefftzWithLf (op, fes, cop_lhs, cop_rhs,
                                                 fes_constraint, nullptr,
                                                 ndof_trefftz, fes_test));
}

/// call `EmbTrefftz` for the plain embedded Trefftz procedure and pack the
/// resulting element matrices in a sparse matrix.
std::tuple<shared_ptr<ngcomp::BaseMatrix>, shared_ptr<ngcomp::BaseVector>>
pythonEmbTrefftzWithLf (shared_ptr<ngfem::SumOfIntegrals> bf,
                        shared_ptr<ngcomp::FESpace> fes,
                        shared_ptr<ngfem::SumOfIntegrals> lf, double eps,
                        shared_ptr<ngcomp::FESpace> test_fes, int tndof,
                        bool getrange, optional<py::dict> stats_dict)
{
  shared_ptr<py::dict> pystats = nullptr;
  if (stats_dict)
    pystats = make_shared<py::dict> (*stats_dict);

  if (fes->IsComplex ())
    {
      std::map<std::string, ngcomp::Vector<Complex>> *stats = nullptr;
      if (pystats)
        stats = new std::map<std::string, ngcomp::Vector<Complex>>;
      auto P = ngcomp::EmbTrefftz<Complex> (bf, fes, lf, eps, test_fes, tndof,
                                            getrange, stats);
      if (pystats)
        for (auto const &x : *stats)
          (*pystats)[py::cast (x.first)] = py::cast (x.second);
      return std::make_tuple (
          ngcomp::Elmats2Sparse<Complex> (std::get<0> (P), fes),
          std::get<1> (P));
    }
  else
    {
      std::map<std::string, ngcomp::Vector<double>> *stats = nullptr;
      if (pystats)
        stats = new std::map<std::string, ngcomp::Vector<double>>;
      auto P = ngcomp::EmbTrefftz<double> (bf, fes, lf, eps, test_fes, tndof,
                                           getrange, stats);
      if (pystats)
        for (auto const &x : *stats)
          (*pystats)[py::cast (x.first)] = py::cast (x.second);
      return std::make_tuple (
          ngcomp::Elmats2Sparse<double> (std::get<0> (P), fes),
          std::get<1> (P));
    }
}

/// call `EmbTrefftz` for the plain embedded Trefftz procedure and pack the
/// resulting element matrices in a sparse matrix.
shared_ptr<ngcomp::BaseMatrix>
pythonEmbTrefftz (shared_ptr<ngfem::SumOfIntegrals> bf,
                  shared_ptr<ngcomp::FESpace> fes, double eps,
                  shared_ptr<ngcomp::FESpace> test_fes, int tndof,
                  bool getrange, optional<py::dict> stats_dict)
{
  shared_ptr<py::dict> pystats = nullptr;
  if (stats_dict)
    pystats = make_shared<py::dict> (*stats_dict);

  if (fes->IsComplex ())
    {
      std::map<std::string, ngcomp::Vector<Complex>> *stats = nullptr;
      if (pystats)
        stats = new std::map<std::string, ngcomp::Vector<Complex>>;
      auto P = std::get<0> (ngcomp::EmbTrefftz<Complex> (
          bf, fes, nullptr, eps, test_fes, tndof, getrange, stats));
      if (pystats)
        for (auto const &x : *stats)
          (*pystats)[py::cast (x.first)] = py::cast (x.second);
      return ngcomp::Elmats2Sparse<Complex> (P, fes);
    }
  else
    {
      std::map<std::string, ngcomp::Vector<double>> *stats = nullptr;
      if (pystats)
        stats = new std::map<std::string, ngcomp::Vector<double>>;
      auto P = std::get<0> (ngcomp::EmbTrefftz<double> (
          bf, fes, nullptr, eps, test_fes, tndof, getrange, stats));
      if (pystats)
        for (auto const &x : *stats)
          (*pystats)[py::cast (x.first)] = py::cast (x.second);
      return ngcomp::Elmats2Sparse<double> (P, fes);
    }
}

void ExportEmbTrefftz (py::module m)
{
  ExportETSpace<ngcomp::L2HighOrderFESpace> (m, "L2EmbTrefftzFESpace");
  // ExportETSpace<ngcomp::VectorL2FESpace,
  // shared_ptr<ngcomp::VectorL2FESpace>> ( m, "VL2EmbTrefftzFESpace");
  ExportETSpace<ngcomp::MonomialFESpace> (m, "MonomialEmbTrefftzFESpace");
  ExportETSpace<ngcomp::CompoundFESpace> (m, "CompoundEmbTrefftzFESpace");

  m.def (
      "EmbeddedTrefftzFES",
      [] (shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::FESpace> {
        shared_ptr<ngcomp::FESpace> nfes;
        if (dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::L2HighOrderFESpace>> (
              dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes));
        // else if (dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes))
        // nfes = make_shared<ngcomp::EmbTrefftzFESpace<
        // ngcomp::VectorL2FESpace, shared_ptr<ngcomp::VectorL2FESpace>>> (
        // dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes));
        else if (dynamic_pointer_cast<ngcomp::MonomialFESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::MonomialFESpace>> (
              dynamic_pointer_cast<ngcomp::MonomialFESpace> (fes));
        else if (dynamic_pointer_cast<ngcomp::CompoundFESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::CompoundFESpace>> (
              dynamic_pointer_cast<ngcomp::CompoundFESpace> (fes));
        else
          throw Exception ("Unknown base fes");
        return nfes;
      },
      R"mydelimiter(
                        Given a FESpace this wrapper produces a Trefftz FESpace using local projections, following the Embedded Trefftz-DG methodology. Use EmbTrefftzFES.SetOp() to set the operator used to construct the embedding.

                        :param fes: FESpace to be wrapped.

                        :return: EmbTrefftzFES
                        )mydelimiter",
      py::arg ("fes"));

  m.def ("TrefftzEmbedding", &pythonEmbTrefftzWithLf,
         R"mydelimiter(
                Computes the Trefftz embedding and particular solution.

                :param bf: operator for which the Trefftz embedding is computed.
                :param fes: DG finite element space of the weak formulation.
                :param lf: Rhs used to compute the particular solution.
                :param eps: Threshold for singular values to be considered zero, defaults to 0
                :param test_fes: Used if test space differs from trial space, defaults to None
                :param tndof: If known, local ndofs of the Trefftz space, else eps and/or test_fes are used to find the dimension
                :param getrange: If True, extract the range instead of the kernel
                :param stats_dict: Pass a dictionary to fill it with stats on the singular values.

                :return: [Trefftz embedding, particular solution]
            )mydelimiter",
         py::arg ("bf"), py::arg ("fes"), py::arg ("lf"), py::arg ("eps") = 0,
         py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0,
         py::arg ("getrange") = false, py::arg ("stats_dict") = nullopt);

  m.def ("TrefftzEmbedding", &pythonEmbTrefftz,
         R"mydelimiter(
                Used without the parameter lf as input the function only returns the Trefftz embedding.

                :return: Trefftz embedding
            )mydelimiter",
         py::arg ("bf"), py::arg ("fes"), py::arg ("eps") = 0,
         py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0,
         py::arg ("getrange") = false, py::arg ("stats_dict") = py::none ());

  m.def ("TrefftzEmbedding", &pythonConstrTrefftz,
         R"mydelimiter(
                creates an embedding matrix P for the given operations `op`,
                `cop_lhs`, `cop_rhs`.
                The embedding is subject to the constraints specified in
                `cop_lhs` and `cop_rhs`.

                 :param op: the differential operation. Can be None
                 :param fes: the finite element space of `op`
                 :param cop_lhs: left hand side of the constraint operation
                 :param cop_rhs: right hand side of the constraint operation
                 :param fes_constraint: finite element space of the constraint operation
                 :param ndof_trefftz: number of degrees of freedom per element
                     in the Trefftz finite element space on `fes`, generated by `op`
                     (i.e. the local dimension of the kernel of `op` on one element)

                 :return: P, the embedding matrix.
   )mydelimiter",
         py::arg ("op"), py::arg ("fes"), py::arg ("cop_lhs"),
         py::arg ("cop_rhs"), py::arg ("fes_constraint"),
         py::arg ("ndof_trefftz"), py::arg ("fes_test") = py::none ());

  m.def ("TrefftzEmbedding", &pythonConstrTrefftzWithLf,
         R"mydelimiter(
                creates an embedding matrix P for the given operations `op`,
                `cop_lhs`, `cop_rhs`.
                The embedding is subject to the constraints specified in
                `cop_lhs` and `cop_rhs`.
                Also generates a particular solution `u_lf`.

                 :param op: the differential operation. Can be None
                 :param fes: the finite element space of `op`
                 :param cop_lhs: left hand side of the constraint operation
                 :param cop_rhs: right hand side of the constraint operation
                 :param fes_constraint: finite element space of the constraint operation
                 :param linear_form: right hand side of the var. formulation
                 :param ndof_trefftz: number of degrees of freedom per element
                     in the Trefftz finite element space on `fes`, generated by `op`
                     (i.e. the local dimension of the kernel of `op` on one element)

                 :return: P, the embedding matrix.
   )mydelimiter",
         py::arg ("op"), py::arg ("fes"), py::arg ("cop_lhs"),
         py::arg ("cop_rhs"), py::arg ("fes_constraint"),
         py::arg ("linear_form"), py::arg ("ndof_trefftz"),
         py::arg ("fes_test") = py::none ());
}

#endif // NGS_PYTHON
