#include "embtrefftz.hpp"
#include "monomialfespace.hpp"
#include <cfloat>
#include <compressedfespace.hpp>

using namespace ngbla;
using namespace ngcomp;

/// @returns `max(a-b, 0)` without integer underflow risk, if `b > a`.
INLINE size_t max_a_minus_b_or_0 (const size_t a, const size_t b) noexcept
{
  return (a >= b) ? a - b : 0;
}

/// Derives the correct local Trefftz ndof form several parameters.
template <typename SCAL>
size_t calcNdofTrefftz (const size_t ndof, const size_t ndof_test,
                        const size_t ndof_conforming,
                        const std::variant<size_t, double> ndof_trefftz,
                        const bool trefftz_op_is_null,
                        const SliceVector<SCAL> singular_values)
{

  if (trefftz_op_is_null)
    {
      // ndof - ndof_conforming might underflow
      return max_a_minus_b_or_0 (ndof, ndof_conforming);
    }
  else if (holds_alternative<size_t> (ndof_trefftz))
    {
      return std::get<size_t> (ndof_trefftz);
    }
  else
    {
      const double eps = std::get<double> (ndof_trefftz);

      // ndof - (ndof_test + ndof_conforming) might underflow
      size_t tndof = max_a_minus_b_or_0 (ndof, ndof_test + ndof_conforming);

      for (const SCAL sigma : singular_values)
        {
          if (abs (sigma) < eps)
            tndof++;
        }
      return tndof;
    }
}

/// @throw std::invalid_argument The number of columns in the matrix must be
///   equal to the number of DofIds in `dof_nrs`
template <typename SCAL, typename TDIST>
void reorderMatrixColumns (
    MatrixView<SCAL, ngbla::RowMajor, size_t, size_t, TDIST> &matrix,
    const Array<DofId> &dof_nrs, LocalHeap &lh)
{
  const auto [n, m] = matrix.Shape ();
  if (m != dof_nrs.Size ())
    throw std::invalid_argument (
        "the width of the matrix must match the length of the dof_nrs");

  const auto heap_reset = HeapReset (lh);

  FlatArray<int> map (m, lh);
  for (size_t i = 0; i < map.Size (); i++)
    map[i] = i;

  QuickSortI (dof_nrs, map);
  auto matrix_copy = FlatMatrix<SCAL> (n, m, lh);
  matrix_copy = matrix;
  for (auto j : Range (m))
    matrix.Col (j) = matrix_copy.Col (map[j]);
}

template <typename SCAL, typename TDIST>
inline void addIntegrationToElementMatrix (
    MatrixView<SCAL, RowMajor, size_t, size_t, TDIST> elmat,
    const Array<shared_ptr<BilinearFormIntegrator>> &bf_integrators,
    const MeshAccess &ma, const ElementId &element_id, const FESpace &fes,
    const FESpace &fes_test, LocalHeap &lh)
{
  const HeapReset hr (lh);

  auto &trafo = ma.GetTrafo (element_id, lh);

  auto &test_fel = fes_test.GetFE (element_id, lh);
  auto &trial_fel = fes.GetFE (element_id, lh);

  const bool mixed_mode = std::addressof (fes_test) != std::addressof (fes);

  bool symmetric_so_far = true;
  size_t bfi_ind = 0;
  while (bfi_ind < bf_integrators.Size ())
    {
      auto &bfi = bf_integrators[bfi_ind];
      bfi_ind++;
      if (bfi->DefinedOnElement (element_id.Nr ()))
        {
          auto &mapped_trafo
              = trafo.AddDeformation (bfi->GetDeformation ().get (), lh);
          try
            {
              if (mixed_mode)
                {
                  const auto &mixed_fel
                      = MixedFiniteElement (trial_fel, test_fel);
                  bfi->CalcElementMatrixAdd (mixed_fel, mapped_trafo, elmat,
                                             symmetric_so_far, lh);
                }
              else
                {
                  bfi->CalcElementMatrixAdd (test_fel, mapped_trafo, elmat,
                                             symmetric_so_far, lh);
                }
            }
          catch (ExceptionNOSIMD const &e)
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

/// Extracts the submatrix of `elmat`,
/// consisting only of entries `elmat[i,j]`, where `i` and `j`
/// are visible dofs.
///
/// Additionally, `dofs` and `dofs_test` are overwritten by only the visible
/// dofs.
///
/// @returns the extracted submatrix, allocated on the `LocalHeap`.
template <typename SCAL, typename TDIST>
FlatMatrix<SCAL> extractVisibleDofs (
    const MatrixView<SCAL, RowMajor, size_t, size_t, TDIST> &elmat,
    const ElementId &element_id, const FESpace &fes, const FESpace &fes_test,
    const shared_ptr<const FESpace> fes_conformity, Array<DofId> &dofs,
    Array<DofId> &dofs_test, Array<DofId> &conformity_dofs, LocalHeap &lh,
    bool compute_new_dofs = false, shared_ptr<BitArray> ignoredofs = nullptr)
{
  Array<DofId> vdofs, vdofs_test, vdofs_conformity;

  fes.GetDofNrs (element_id, vdofs, VISIBLE_DOF);
  if (ignoredofs)
    for (size_t i = 0; i < vdofs.Size (); i++)
      if (ignoredofs->Test (vdofs[i]))
        vdofs.RemoveElement (i);
  // vdofs.DeleteElement (i); if we dont need to retain the order
  fes_test.GetDofNrs (element_id, vdofs_test, VISIBLE_DOF);
  if (fes_conformity)
    fes_conformity->GetDofNrs (element_id, vdofs_conformity, VISIBLE_DOF);

  FlatMatrix<SCAL> velmat (vdofs_test.Size () + vdofs_conformity.Size (),
                           vdofs.Size (), lh);

  size_t conformity_offset = fes_conformity ? conformity_dofs.Size () : 0;
  size_t vconformity_offset = fes_conformity ? vdofs_conformity.Size () : 0;

  if (fes_conformity)
    {
      for (size_t vi = 0; vi < vdofs_conformity.Size (); vi++)
        {
          const size_t i = conformity_dofs.Pos (vdofs_conformity[vi]);
          for (size_t vj = 0; vj < vdofs.Size (); vj++)
            {
              const size_t j = dofs.Pos (vdofs[vj]);
              velmat (vi, vj) = elmat (i, j);
            }
        }
    }
  for (size_t vi = 0; vi < vdofs_test.Size (); vi++)
    {
      const size_t i = dofs_test.Pos (vdofs_test[vi]);
      for (size_t vj = 0; vj < vdofs.Size (); vj++)
        {
          const size_t j = dofs.Pos (vdofs[vj]);
          velmat (vconformity_offset + vi, vj)
              = elmat (conformity_offset + i, j);
        }
    }

  if (compute_new_dofs)
    {
      conformity_dofs = std::move (vdofs_conformity);
      dofs = std::move (vdofs);
      dofs_test = std::move (vdofs_test);
    }
  return velmat;
}

template <typename SCAL, typename TDIST>
Matrix<SCAL> putbackVisibleDofs (
    const MatrixView<SCAL, RowMajor, size_t, size_t, TDIST> &velmat,
    const ElementId &element_id, const FESpace &fes, const Array<DofId> &vdofs,
    shared_ptr<BitArray> ignoredofs = nullptr)
{
  Array<DofId> dofs;
  fes.GetDofNrs (element_id, dofs);

  size_t all_ignored_dofs = 0;
  if (ignoredofs)
    for (size_t i = 0; i < dofs.Size (); i++)
      if (ignoredofs->Test (dofs[i]))
        all_ignored_dofs++;

  Matrix<SCAL> elmat (dofs.Size (), velmat.Width () + all_ignored_dofs);
  elmat = static_cast<SCAL> (0.0);

  for (size_t j = 0, ignored_dofs = 0; j < dofs.Size (); j++)
    {
      const size_t vj = vdofs.Pos (dofs[j]);
      if (vj != size_t (-1))
        // elmat.Row (j) = velmat.Row (vj);
        for (size_t i = 0; i < velmat.Width (); i++)
          elmat (j, all_ignored_dofs + i) = velmat (vj, i);
      if (ignoredofs && ignoredofs->Test (dofs[j]))
        {
          elmat (j, ignored_dofs++) = 1.0;
        }
    }
  return elmat;
}

void fesFromOp (const SumOfIntegrals &op, shared_ptr<FESpace> &fes,
                shared_ptr<FESpace> &fes_test)
{
  Array<shared_ptr<Integral>> ics;

  for (shared_ptr<Integral> icf : op.icfs)
    icf->cf->TraverseTree ([&] (CoefficientFunction &nodecf) {
      auto proxy = dynamic_pointer_cast<ProxyFunction> (
          (&nodecf)->shared_from_this ());
      if (proxy && proxy->IsTrialFunction () && !fes)
        fes = proxy->GetFESpace ();
      else if (proxy && proxy->IsTrialFunction ()
               && fes != proxy->GetFESpace ())
        throw Exception ("Two different trial FESpaces in the same operator");
      else if (proxy && proxy->IsTestFunction () && !fes_test)
        fes_test = proxy->GetFESpace ();
      else if (proxy && proxy->IsTestFunction ()
               && fes_test != proxy->GetFESpace ())
        throw Exception ("Two different test FESpaces in the same operator");
    });
}

/// Fills the two creators with the sparsity pattern needed for
/// the Trefftz embedding.
/// @tparam NZ_FUNC has signature `(ElementId) -> size_t`
template <typename SCAL>
size_t
createConformingTrefftzTables (Table<int> &table, Table<int> &table2,
                               const FlatArray<optional<Matrix<SCAL>>> &etmats,
                               const FlatArray<size_t> &local_ndofs_trefftz,
                               const FESpace &fes,
                               shared_ptr<const FESpace> fes_conformity,
                               shared_ptr<const BitArray> ignoredofs)
{
  const auto ma = fes.GetMeshAccess ();
  const size_t ne = ma->GetNE (VOL);
  const size_t ndof_conforming
      = (fes_conformity) ? fes_conformity->GetNDof () : 0;
  TableCreator<int> creator (ne);
  TableCreator<int> creator2 (ne);

  size_t trefftz_dof_offset = ndof_conforming;

  Array<DofId> new_ignore_dofnrs; // TODO: should be keyd
  if (ignoredofs)
    {
      new_ignore_dofnrs.SetSize (fes.GetNDof ());
      for (size_t i = 0; i < fes.GetNDof (); i++)
        if (ignoredofs->Test (i))
          new_ignore_dofnrs[i] = trefftz_dof_offset++;
    }

  size_t global_trefftz_ndof = 0;
  for (; !creator.Done (); creator++, creator2++)
    {
      size_t next_trefftz_dof = trefftz_dof_offset;
      for (auto ei : ma->Elements (VOL))
        {
          if (!etmats[ei.Nr ()])
            continue;

          // first compute the Trefftz dofs. The dof numbers of the Trefftz
          // dofs are shifted up by conforming_ndof, to avoid conflicts between
          // Trefftz and Constraint dofs.
          size_t nz = local_ndofs_trefftz[ei.Nr ()];
          Array<DofId> dnums;
          fes.GetDofNrs (ei, dnums);
          bool hasregdof = false;
          for (DofId d : dnums)
            if (IsRegularDof (d))
              {
                creator.Add (ei.Nr (), d);
                hasregdof = true;
                if (ignoredofs && ignoredofs->Test (d))
                  creator2.Add (ei.Nr (), new_ignore_dofnrs[d]);
              }
          // assumption here: Either all or no dof is regular
          if (hasregdof)
            {
              for (size_t d = 0; d < nz; d++)
                creator2.Add (ei.Nr (), next_trefftz_dof++);
            }

          // then compute the Constraint dofs.
          if (fes_conformity)
            {
              Array<DofId> dofs_conforming;
              fes_conformity->GetDofNrs (ei, dofs_conforming, VISIBLE_DOF);

              bool hasregdof = false;
              for (DofId d : dofs_conforming)
                if (IsRegularDof (d))
                  {
                    // creator.Add (ei.Nr (), d);
                    hasregdof = true;
                  }
              // assumption here: Either all or no dof is regular
              if (hasregdof)
                {
                  for (DofId d : dofs_conforming)
                    creator2.Add (ei.Nr (), d);
                }
            }
        }
      global_trefftz_ndof = next_trefftz_dof - ndof_conforming;
    }

  (*testout) << "created " << global_trefftz_ndof << " many trefftz dofs"
             << std::endl;

  table = creator.MoveTable ();
  table2 = creator2.MoveTable ();
  return global_trefftz_ndof + ndof_conforming;
}

template <typename SCAL>
INLINE void
fillSparseMatrixWithData (SparseMatrix<SCAL> &P,
                          const FlatArray<optional<Matrix<SCAL>>> &etmats,
                          const Table<int> &table, const Table<int> &table2,
                          const MeshAccess &ma,
                          shared_ptr<BitArray> ignoredofs)
{
  P.SetZero ();
  for (auto ei : ma.Elements (VOL))
    if (etmats[ei.Nr ()])
      {
        P.AddElementMatrix (table[ei.Nr ()], table2[ei.Nr ()],
                            *etmats[ei.Nr ()]);
      }

  if (ignoredofs)
    {
      Array<int> id_index (ignoredofs->NumSet ());
      for (size_t i = 0, id = 0; i < ignoredofs->Size (); i++)
        if (ignoredofs->Test (i))
          {
            id_index[id] = P.GetPosition (i, id);
            id++;
          }
      Vector<SCAL> one (id_index.Size ());
      one = 1.0;
      P.AsVector ().SetIndirect (id_index, one);
    }
}

namespace ngcomp
{
  /// assembles a global sparse matrix from the given element matrices.
  /// @param etmats vector of all element matrices
  /// @param fes non-Trefftz finite element space
  /// @tparam SCAL scalar type of the matrix entries
  template <typename SCAL>
  shared_ptr<BaseMatrix>
  Elmats2Sparse (const FlatArray<optional<Matrix<SCAL>>> &etmats,
                 const FlatArray<size_t> &local_ndofs_trefftz,
                 const FESpace &fes, shared_ptr<const FESpace> fes_conformity,
                 shared_ptr<BitArray> ignoredofs)
  {
    const auto ma = fes.GetMeshAccess ();

    Table<int> table, table2;
    const size_t conformity_plus_trefftz_dofs = createConformingTrefftzTables (
        table, table2, etmats, local_ndofs_trefftz, fes, fes_conformity,
        ignoredofs);

    auto P = make_shared<SparseMatrix<SCAL>> (
        fes.GetNDof (), conformity_plus_trefftz_dofs, table, table2, false);
    fillSparseMatrixWithData (*P, etmats, table, table2, *ma, ignoredofs);

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
    const SumOfIntegrals &trhs,
    Array<shared_ptr<LinearFormIntegrator>> lfis[4])
{
  for (auto icf : trhs.icfs)
    {
      DifferentialSymbol &dx = icf->dx;
      lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
    }
}

/// @returns true, if `fes` is a CompressedFESpace
/// or is a CompoundFESpace with a CompressedFESpace as a child.
bool fesHasCompressedChild (const FESpace &fes)
{
  if (dynamic_cast<const CompressedFESpace *> (&fes))
    return true;

  if (dynamic_cast<const CompoundFESpace *> (&fes))
    {
      for (const shared_ptr<FESpace> &component_fes :
           dynamic_cast<const CompoundFESpace &> (fes).Spaces ())
        {
          if (fesHasCompressedChild (*component_fes))
            return true;
        }
    }
  return false;
}

/// Determines, if the FESpace has unused or hidden dofs
/// (which are considered inactive for this purpose).
///
/// If the provided space is a CompressedFESpace,
/// it is assumed that the base space will have inactive dofs.
bool fesHasInactiveDofs (const FESpace &fes)
{
  if (fesHasCompressedChild (fes))
    return true;
  const size_t ndof = fes.GetNDof ();
  for (size_t d = 0; d < ndof; d++)
    if (fes.GetDofCouplingType (d) <= HIDDEN_DOF)
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

  /// `A = U * Sigma * V`
  /// A gets overwritten with Sigma
  /// @param A has dimension n * m
  /// @param U has dimension n * n
  /// @param V has dimension m * m
  template <typename SCAL, typename TDIST>
  void getSVD (MatrixView<SCAL, ngbla::RowMajor, size_t, size_t, TDIST> A,
               FlatMatrix<SCAL, ColMajor> U, FlatMatrix<SCAL, ColMajor> V)
  {
    auto [height, width] = A.Shape ();
    Matrix<SCAL, ColMajor> AA (height, width);
    AA = A;

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
    for (size_t i = 0; i < min (A.Width (), A.Height ()); i++)
      A (i, i) = AA (i, i);
  }
}

/// @param num_zero number of singular values `sigma_i = 0` in the matrix Sigma
/// @return pseudoinverse of A (as some `ngbla::Expr` type to avoid
/// allocations)
template <typename SCAL, typename TDIST>
Matrix<SCAL>
invertSVD (const FlatMatrix<SCAL, ColMajor> &U,
           const MatrixView<SCAL, RowMajor, size_t, size_t, TDIST> &Sigma,
           const FlatMatrix<SCAL, ColMajor> &V, size_t num_zero, LocalHeap &lh)
{
  const HeapReset hr (lh);
  const auto [m, n] = Sigma.Shape ();

  FlatMatrix<SCAL> sigma_inv_times_ut (n, U.Height (), lh);
  for (size_t i = 0; i < n; i++)
    if (i < min (n - num_zero, m))
      sigma_inv_times_ut.Row (i) = 1.0 / Sigma (i, i) * U.Col (i);
    else
      sigma_inv_times_ut.Row (i) = SCAL (0.);

  return Trans (V) * sigma_inv_times_ut;
}

/// calculates from the given space and linear form integrators a particular
/// solution. note: allocates the solution on the given local heap.
/// @param lfis arrays of linear form integrators, for `VOL`, `BND`, `BBND`,
/// `BBBND`
/// @return the element-local solution, allocated on the local heap `mlh`
/// @tparam T type of matrix expression
template <typename SCAL, typename T>
void calculateParticularSolution (
    FlatVector<SCAL> partsol, Array<shared_ptr<LinearFormIntegrator>> lfis[4],
    const FESpace &fes_test, const ElementId ei, const MeshAccess &ma,
    const Expr<T> &inverse_elmat, LocalHeap &mlh)
{
  Array<DofId> dofs_test;
  fes_test.GetDofNrs (ei, dofs_test);

  const HeapReset hr (mlh);
  const auto &test_fel = fes_test.GetFE (ei, mlh);
  const auto &trafo = ma.GetTrafo (ei, mlh);
  FlatVector<SCAL> elvec (dofs_test.Size (), mlh),
      elveci (dofs_test.Size (), mlh);
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

  if (fesHasInactiveDofs (fes_test))
    {
      Array<DofId> vdofs_test;
      fes_test.GetDofNrs (ei, vdofs_test, VISIBLE_DOF);
      FlatVector<SCAL> velvec (vdofs_test.Size (), mlh);

      for (size_t vi = 0; vi < vdofs_test.Size (); vi++)
        {
          const size_t i = dofs_test.Pos (vdofs_test[vi]);
          velvec[vi] = elvec[i];
        }
      // elvec = std::move(velvec);
      elvec.AssignMemory (velvec.Size (), &velvec[0]);
    }
  partsol = inverse_elmat * elvec;
}

namespace ngcomp
{
  mutex stats_mutex;

  template <typename SCAL> void TrefftzEmbedding::EmbTrefftz ()
  {
    static_assert (std::is_same_v<double, SCAL>
                       || std::is_same_v<Complex, SCAL>,
                   "SCAL must be double or Complex");

    // statistics stuff
    Vector<double> sing_val_avg;
    Vector<double> sing_val_max;
    Vector<double> sing_val_min;
    std::atomic<size_t> active_elements = 0;

    auto ma = fes->GetMeshAccess ();
    const size_t num_elements = ma->GetNE (VOL);
    // #TODO what is a good size for the local heap?
    // For the moment: large enough constant size.
    LocalHeap clh = LocalHeap (100 * 1000 * 1000, "embt", true);

    // calculate the integrators for the three bilinear forms,
    // each for VOL, BND, BBND, BBBND, hence 4 arrays per bilnear form
    Array<shared_ptr<BilinearFormIntegrator>> op_integrators[4],
        cop_lhs_integrators[4], cop_rhs_integrators[4];
    if (top)
      calculateBilinearFormIntegrators (*top, op_integrators);
    if (cop && crhs)
      {
        calculateBilinearFormIntegrators (*cop, cop_lhs_integrators);
        calculateBilinearFormIntegrators (*crhs, cop_rhs_integrators);
      }

    Array<optional<Matrix<SCAL>>> etmats (num_elements);
    Array<optional<Matrix<SCAL>>> etmats_inv (num_elements);
    Array<size_t> local_ndofs_trefftz (num_elements);

    const bool any_fes_has_inactive_dofs
        = fesHasInactiveDofs (*fes) || fesHasInactiveDofs (*fes_test)
          || ((fes_conformity) ? fesHasInactiveDofs (*fes_conformity) : false);

    auto particular_solution_vec
        = make_shared<VVector<SCAL>> (fes->GetNDof ());
    particular_solution_vec->operator= (0.0);
    Array<shared_ptr<LinearFormIntegrator>> lfis[4];
    if (trhs)
      calculateLinearFormIntegrators (*trhs, lfis);

    // we compute embedding matrices E_C, E_T and particular solution u_p that
    // satisfy
    //  /C_l\  * (E_C | E_T | u_p) = / C_r | 0 | 0 \                       //
    //  \ L /                        \  0  | 0 | f /                       //
    // C_l.shape == (ndof_conforming, ndof),
    // L.shape == (ndof_test, ndof)
    // C_r.shape == (ndof_conforming, ndof_conforming),
    // E_C.shape == (ndof, ndof_conforming),
    // E_T.shape == (ndof, ndof_trefftz),
    // f.shape == (ndof_test), u_p.shape == (ndof)
    ma->IterateElements (
        VOL, clh, [&] (Ngs_Element mesh_element, LocalHeap &lh) {
          const ElementId element_id = ElementId (mesh_element);

          // skip this element, if the bilinear forms are not defined
          // on this element
          if (!(top && bfIsDefinedOnElement (*top, mesh_element))
              && !(cop && bfIsDefinedOnElement (*cop, mesh_element))
              && !(crhs && bfIsDefinedOnElement (*crhs, mesh_element)))
            return;

          Array<DofId> dofs, dofs_test, dofs_conforming;
          fes->GetDofNrs (element_id, dofs);
          fes_test->GetDofNrs (element_id, dofs_test);
          if (fes_conformity)
            fes_conformity->GetDofNrs (element_id, dofs_conforming);

          size_t ndof = dofs.Size ();
          size_t ndof_test = dofs_test.Size ();
          size_t ndof_conforming = dofs_conforming.Size ();
          FlatMatrix<SCAL> elmat_A (ndof_test + ndof_conforming, ndof, lh);
          auto [elmat_Cl, elmat_L] = elmat_A.SplitRows (ndof_conforming);

          FlatMatrix<SCAL> elmat_B (ndof_test + ndof_conforming,
                                    ndof_conforming, lh);
          elmat_A = static_cast<SCAL> (0.);
          elmat_B = static_cast<SCAL> (0.);

          // elmat_cr is a view into elamt_b
          auto elmat_Cr = elmat_B.Rows (ndof_conforming);

          // the diff. operator L operates only on volume terms
          addIntegrationToElementMatrix (elmat_L, op_integrators[VOL], *ma,
                                         element_id, *fes, *fes_test, lh);
          if (fes_conformity)
            {
              for (const auto vorb : { VOL, BND, BBND, BBBND })
                {
                  addIntegrationToElementMatrix (
                      elmat_Cl, cop_lhs_integrators[vorb], *ma, element_id,
                      *fes, *fes_conformity, lh);
                  addIntegrationToElementMatrix (
                      elmat_Cr, cop_rhs_integrators[vorb], *ma, element_id,
                      *fes_conformity, *fes_conformity, lh);
                }
            }

          if (any_fes_has_inactive_dofs || ignoredofs)
            {
              if (fes_conformity)
                {
                  elmat_B.Assign (extractVisibleDofs (
                      elmat_B, element_id, *fes_conformity, *fes_test,
                      fes_conformity, dofs_conforming, dofs_test,
                      dofs_conforming, lh));
                }

              elmat_A.Assign (extractVisibleDofs (
                  elmat_A, element_id, *fes, *fes_test, fes_conformity, dofs,
                  dofs_test, dofs_conforming, lh, true, ignoredofs));

              ndof = dofs.Size ();
              ndof_test = dofs_test.Size ();
              ndof_conforming = dofs_conforming.Size ();

              // Technically not needed anymore,
              // but since elmat_Cl, elmat_L, and elmat_Cr are still living
              // (but now pointing to the wrong memory), this is an
              // easy-to-miss source of errors in the future.
              auto [elmat_Cl_tmp, elmat_L_tmp]
                  = elmat_A.SplitRows (ndof_conforming);
              elmat_Cl.Assign (elmat_Cl_tmp);
              elmat_L.Assign (elmat_L_tmp);
              elmat_Cr.Assign (elmat_B.Rows (ndof_conforming));
            }
          // reorder elmat_cr. Needs to happen after visible dof extraction.
          // #TODO: why is this necessary?
          reorderMatrixColumns (elmat_Cr, dofs_conforming, lh);

          FlatMatrix<SCAL, ColMajor> U (elmat_A.Height (), lh),
              V (elmat_A.Width (), lh);
          getSVD<SCAL> (elmat_A, U, V);

          // # TODO: incorporate the double variant
          const size_t ndof_trefftz_i
              = calcNdofTrefftz (ndof, ndof_test, ndof_conforming,
                                 ndof_trefftz, !top, elmat_A.Diag (0));
          if (ndof_trefftz_i + ndof_conforming == 0)
            throw std::invalid_argument ("zero trefftz dofs");

          const auto elmat_A_inv_expr
              = invertSVD (U, elmat_A, V, ndof_trefftz_i, lh);
          FlatMatrix<SCAL> elmat_A_inv (ndof, ndof_conforming + ndof_test, lh);
          // Calculate the matrix entries and write them to memory.
          elmat_A_inv = elmat_A_inv_expr;

          // T = (T_c | T_t)
          Matrix<SCAL> elmat_T (ndof, ndof_trefftz_i + ndof_conforming);
          auto [elmat_Tc, elmat_Tt] = elmat_T.SplitCols (ndof_conforming);

          // T_c solves A @ T_c = B,
          elmat_Tc = elmat_A_inv * elmat_B;

          // if (get_range)
          // elmat_Tt = U.Cols (0, dofs.Size () - ndof_trefftz_i);
          elmat_Tt = Trans (V.Rows (ndof - ndof_trefftz_i, ndof));

          if (trhs)
            {
              auto elmat_T_inv = elmat_A_inv.Cols (
                  ndof_conforming, ndof_conforming + ndof_test);
              FlatVector<SCAL> partsol (dofs.Size (), lh);
              calculateParticularSolution<SCAL> (
                  partsol, lfis, *fes_test, element_id, *ma, elmat_T_inv, lh);
              particular_solution_vec->SetIndirect (dofs, partsol);
            }

          if (ignoredofs)
            elmat_T = putbackVisibleDofs (elmat_T, element_id, *fes, dofs,
                                          ignoredofs);

          etmats[element_id.Nr ()] = make_optional<Matrix<SCAL>> (elmat_T);
          etmats_inv[element_id.Nr ()]
              = make_optional<Matrix<SCAL>> (elmat_A_inv);
          local_ndofs_trefftz[element_id.Nr ()] = ndof_trefftz_i;

          if (stats)
            {
              const lock_guard<mutex> lock (stats_mutex);
              if (sing_val_avg.Size () == 0)
                {
                  sing_val_avg.SetSize (elmat_A.Height ());
                  sing_val_max.SetSize (elmat_A.Height ());
                  sing_val_min.SetSize (elmat_A.Height ());
                  sing_val_avg = 0;
                  sing_val_max = 0;
                  sing_val_min = DBL_MAX;
                }
              active_elements += 1;
              for (size_t i = 0; i < elmat_A.Height (); i++)
                {
                  sing_val_avg[i] += abs (elmat_A (i, i));
                  sing_val_max[i]
                      = max (sing_val_max[i], abs (elmat_A (i, i)));
                  sing_val_min[i]
                      = min (sing_val_min[i], abs (elmat_A (i, i)));
                }
            }

          (*testout) << "element " << element_id << endl
                     << "fes has ndof:" << ndof
                     << "fes_test has ndof:" << ndof_test
                     << "fes_conformity has ndof:" << ndof_conforming << endl
                     << "elmat_t1" << endl
                     << elmat_Tc << endl
                     << "elmat_t2" << endl
                     << elmat_Tt << endl;
        });
    if (stats)
      {
        sing_val_avg /= active_elements.load ();
        (*stats)["singavg"] = Vector<double> (sing_val_avg);
        (*stats)["singmax"] = Vector<double> (sing_val_max);
        (*stats)["singmin"] = Vector<double> (sing_val_min);
      }

    if constexpr (std::is_same_v<double, SCAL>)
      {
        this->etmats = etmats;
        this->etmats_inv = etmats_inv;
      }
    else
      {
        this->etmatsc = etmats;
        this->etmatsc_inv = etmats_inv;
      }
    this->local_ndofs_trefftz = local_ndofs_trefftz;
    this->psol = particular_solution_vec;
  }

  TrefftzEmbedding::TrefftzEmbedding (
      shared_ptr<SumOfIntegrals> _top, shared_ptr<SumOfIntegrals> _trhs,
      shared_ptr<SumOfIntegrals> _cop, shared_ptr<SumOfIntegrals> _crhs,
      size_t _ndof_trefftz, double _eps, shared_ptr<FESpace> _fes,
      shared_ptr<FESpace> _fes_test, shared_ptr<FESpace> _fes_conformity,
      shared_ptr<BitArray> _ignoredofs)
      : top (_top), trhs (_trhs), cop (_cop), crhs (_crhs),
        ignoredofs (_ignoredofs)
  {
    if (_ndof_trefftz == 0)
      ndof_trefftz = _eps;
    else if (_eps == 0.0)
      ndof_trefftz = _ndof_trefftz;
    else
      throw std::invalid_argument (
          "cannot use ndof_trefftz and eps at the same time");

    if (top)
      fesFromOp (*top, fes, fes_test);
    if (cop)
      fesFromOp (*cop, fes, fes_conformity);
    if (_fes)
      fes = _fes;
    if (_fes_test)
      fes_test = _fes_test;
    if (_fes_conformity)
      fes_conformity = _fes_conformity;
    if (!fes_test)
      fes_test = fes;
    if (!fes)
      throw Exception ("TrefftzEmbedding: no fespace found");

    ma = fes->GetMeshAccess ();

    if (fes->IsComplex ())
      EmbTrefftz<Complex> ();
    else
      EmbTrefftz<double> ();
  }

  shared_ptr<const BaseVector> TrefftzEmbedding::GetParticularSolution () const
  {
    return psol;
  }

  shared_ptr<const BaseVector> TrefftzEmbedding::GetParticularSolution (
      shared_ptr<SumOfIntegrals> _trhs) const
  {
    LocalHeap clh = LocalHeap (100 * 1000 * 1000, "embt_psol", true);

    shared_ptr<BaseVector> particular_solution_vec;
    if (fes->IsComplex ())
      particular_solution_vec
          = make_shared<VVector<Complex>> (fes->GetNDof ());
    else
      particular_solution_vec = make_shared<VVector<double>> (fes->GetNDof ());
    particular_solution_vec->operator= (0.0);
    Array<shared_ptr<LinearFormIntegrator>> lfis[4];
    calculateLinearFormIntegrators (*_trhs, lfis);

    ma->IterateElements (
        VOL, clh, [&] (Ngs_Element mesh_element, LocalHeap &lh) {
          if (fes->IsComplex () && !etmatsc[mesh_element.Nr ()])
            return;
          if (!fes->IsComplex () && !etmats[mesh_element.Nr ()])
            return;
          const ElementId element_id = ElementId (mesh_element);

          Array<DofId> dofs, dofs_test;
          fes->GetDofNrs (element_id, dofs, VISIBLE_DOF);
          if (ignoredofs)
            for (size_t i = 0; i < dofs.Size (); i++)
              if (ignoredofs->Test (dofs[i]))
                dofs.RemoveElement (i);
          fes_test->GetDofNrs (element_id, dofs_test, VISIBLE_DOF);

          if (fes->IsComplex ())
            {
              const size_t width = (*etmatsc_inv[element_id.Nr ()]).Width ();
              auto elmat_T_inv = (*etmatsc_inv[element_id.Nr ()])
                                     .Cols (width - dofs_test.Size (), width);
              FlatVector<Complex> partsol (dofs.Size (), lh);
              calculateParticularSolution<Complex> (
                  partsol, lfis, *fes_test, element_id, *ma, elmat_T_inv, lh);
              particular_solution_vec->SetIndirect (dofs, partsol);
            }
          else
            {
              const size_t width = (*etmats_inv[element_id.Nr ()]).Width ();
              auto elmat_T_inv = (*etmats_inv[element_id.Nr ()])
                                     .Cols (width - dofs_test.Size (), width);
              FlatVector<double> partsol (dofs.Size (), lh);
              calculateParticularSolution<double> (
                  partsol, lfis, *fes_test, element_id, *ma, elmat_T_inv, lh);
              particular_solution_vec->SetIndirect (dofs, partsol);
            }
        });
    return particular_solution_vec;
  }

  shared_ptr<const BaseMatrix> TrefftzEmbedding::GetEmbedding () const
  {
    if (fes->IsComplex ())
      return Elmats2Sparse<Complex> (etmatsc, local_ndofs_trefftz, *fes,
                                     fes_conformity, ignoredofs);
    else
      return Elmats2Sparse<double> (etmats, local_ndofs_trefftz, *fes,
                                    fes_conformity, ignoredofs);
  }

  shared_ptr<BaseVector>
  TrefftzEmbedding::Embed (const shared_ptr<const BaseVector> tvec) const
  {
    LocalHeap lh (100 * 1000 * 1000, "embt", true);

    shared_ptr<BaseVector> vec;

    Table<int> table, table2;
    if (fes->IsComplex ())
      {
        vec = make_shared<VVector<Complex>> (fes->GetNDof ());
        createConformingTrefftzTables (table, table2, etmatsc,
                                       local_ndofs_trefftz, *fes,
                                       fes_conformity, ignoredofs);
      }
    else
      {
        vec = make_shared<VVector<double>> (fes->GetNDof ());
        createConformingTrefftzTables (table, table2, etmats,
                                       local_ndofs_trefftz, *fes,
                                       fes_conformity, ignoredofs);
      }
    if (fes_conformity)
      vec->SetScalar (0.0);

    ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
      Array<DofId> dofs, tdofs;
      dofs = table[ei.Nr ()];
      tdofs = table2[ei.Nr ()];

      if (fes->IsComplex ())
        {
          FlatVector<Complex> telvec (tdofs.Size (), mlh);
          tvec->GetIndirect (tdofs, telvec);
          FlatVector<Complex> elvec (dofs.Size (), mlh);
          elvec = *etmatsc[ei.Nr ()] * telvec;
          if (fes_conformity)
            vec->AddIndirect (dofs, elvec, true);
          else
            vec->SetIndirect (dofs, elvec);
        }
      else
        {
          FlatVector<> telvec (tdofs.Size (), mlh);
          tvec->GetIndirect (tdofs, telvec);
          FlatVector<> elvec (dofs.Size (), mlh);
          elvec = *etmats[ei.Nr ()] * telvec;
          if (fes_conformity)
            vec->AddIndirect (dofs, elvec, true);
          else
            vec->SetIndirect (dofs, elvec);
        }
    });
    return vec;
  }

  ////////////////////////// EmbTrefftzFESpace ///////////////////////////

  /// Copies `source` to the beginning of `target`.
  ///
  /// For `source.Size() > target.Size()`,
  /// the behaviour is undefined.
  void copyBitArray (const shared_ptr<BitArray> target,
                     const shared_ptr<const BitArray> source)
  {
    assert (source->Size () <= target->Size ()
            && "The target must not be smaller than the source BitArray.");

    assert (source != nullptr
            && "The source BitArray pointer may not be null");
    assert (target != nullptr
            && "The target BitArray pointer may not be null");

    for (size_t i = 0; i < source->Size (); i++)
      {
        if (source->Test (i))
          target->SetBit (i);
        else
          target->Clear (i);
      }
  }

  template <typename T> void EmbTrefftzFESpace<T>::adjustDofsAfterSetOp ()
  {
    static_assert (std::is_base_of_v<FESpace, T>, "T must be a FESpace");

    // #TODO: what about hidden dofs?
    if (!this->IsComplex ())
      {
        Table<DofId> _table_dummy{};
        createConformingTrefftzTables (_table_dummy, elnr_to_dofs, etmats,
                                       emb->GetLocalNodfsTrefftz (), *fes,
                                       fes_conformity, ignoredofs);
        for (size_t i = 0; i < etmats.Size (); i++)
          if (etmats[i])
            QuickSort (elnr_to_dofs[i]);
      }
    else
      {
        Table<DofId> _table_dummy{};
        createConformingTrefftzTables (_table_dummy, elnr_to_dofs, etmatsc,
                                       emb->GetLocalNodfsTrefftz (), *fes,
                                       fes_conformity, ignoredofs);
        for (size_t i = 0; i < etmatsc.Size (); i++)
          if (etmatsc[i])
            QuickSort (elnr_to_dofs[i]);
      }

    const size_t ndof_conformity
        = (fes_conformity) ? fes_conformity->GetNDof () : 0;

    T::Update ();

    size_t ndof_trefftz = 0;
    for (auto ei : this->ma->Elements (VOL))
      {
        // skip this element, if there is no element matrix defined
        if ((this->IsComplex () && !etmatsc[ei.Nr ()])
            || (!this->IsComplex () && !etmats[ei.Nr ()]))
          continue;

        const size_t ndof_trefftz_local
            = emb->GetLocalNodfsTrefftz ()[ei.Nr ()];
        ndof_trefftz += ndof_trefftz_local;
      }

    size_t ignored_dofs = 0;
    if (ignoredofs)
      ignored_dofs = ignoredofs->NumSet ();

    // The conformity dofs might overlap.
    // Overall, they add up to exactly the number of
    // dofs in the conformity space.
    const size_t new_ndof = ignored_dofs + ndof_conformity + ndof_trefftz;
    this->SetNDof (new_ndof);

    this->ctofdof.SetSize (new_ndof);

    // We start the numbering of the dofs with the conformity dofs,
    // and continue with the Trefftz dofs.
    if (ignoredofs)
      for (size_t i = 0, idof = 0; i < ignoredofs->Size (); i++)
        if (ignoredofs->Test (i))
          this->ctofdof[idof++] = fes->GetDofCouplingType (i);
    for (size_t i = ignored_dofs; i < ndof_conformity; i++)
      this->ctofdof[i] = fes_conformity->GetDofCouplingType (i);
    for (size_t i = ignored_dofs + ndof_conformity; i < new_ndof; i++)
      this->ctofdof[i] = LOCAL_DOF;

    T::FinalizeUpdate ();

    // needs previous FinalizeUpdate to construct free_dofs for `this`
    if (fes_conformity)
      copyBitArray (this->GetFreeDofs (), fes_conformity->GetFreeDofs ());
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    // TODO: ignore dofs for BND, BBND, BBBND?
    if (!T::DefinedOn (ei) || ei.VB () != VOL)
      return;
    // 1. Provide the dof nrs of the conforming Trefftz space, that are
    // associated to the element ei.
    const FlatArray<DofId> tdofnrs = elnr_to_dofs[ei.Nr ()];

    // 2. In order to properly hook into ngsolve's assembly routine,
    // we need to provide as many dofs as the underlying space T has.
    // So, we first provide tdofnrs, which may contain fewer dofs than
    // required. The rest of the array dnums, we fill up with the non-regular
    // DofId NO_DOF_NR_CONDENSE, which marks these excess dofs to be condensed
    // out.
    T::GetDofNrs (ei, dnums);
    for (size_t i = 0; i < dnums.Size (); i++)
      if (IsRegularDof (dnums[i]))
        {
          if (i < tdofnrs.Size ())
            dnums[i] = tdofnrs[i];
          else
            dnums[i] = NO_DOF_NR_CONDENSE;
        }
  }

  /// double/complex generic implementation for the methods VTransformM(R | C)
  /// for the embedded Trefftz FESpace.
  template <typename SCAL>
  void etFesVTransformM (SliceMatrix<SCAL> mat, const TRANSFORM_TYPE type,
                         const Matrix<SCAL> elmat)
  {
    static Timer timer ("EmbTrefftz: MTransform");
    RegionTimer reg (timer);

    Matrix<SCAL> temp_mat (mat.Height (), mat.Width ());

    const size_t tndof = elmat.Width ();

    if (type == TRANSFORM_MAT_LEFT)
      {
        temp_mat.Rows (0, tndof) = Trans (elmat) * mat;
        mat = temp_mat;
      }
    else if (type == TRANSFORM_MAT_RIGHT)
      {
        temp_mat.Cols (0, tndof) = mat * elmat;
        mat = temp_mat;
      }
    else if (type == TRANSFORM_MAT_LEFT_RIGHT)
      {
        auto mat_times_elmat = temp_mat.Cols (0, tndof);
        mat_times_elmat = mat * elmat;

        auto mat_upleft = mat.Cols (0, tndof).Rows (0, tndof);
        mat_upleft = Trans (elmat) * mat_times_elmat;
      }
    else
      {
        stringstream err;
        err << "VTransformM is not implemented for TRANSFORM_TYPE " << type;
        throw std::invalid_argument (err.str ());
      }
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformMR (ElementId ei, SliceMatrix<double> mat,
                                      TRANSFORM_TYPE type) const
  {
    const auto elmat = *etmats[ei.Nr ()];
    etFesVTransformM (mat, type, elmat);
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformMC (ElementId ei, SliceMatrix<Complex> mat,
                                      TRANSFORM_TYPE type) const
  {
    const auto elmat = *etmatsc[ei.Nr ()];
    etFesVTransformM (mat, type, elmat);
  }

  /// double/complex generic implementation for the methods VTransformV(R | C)
  /// for the embedded Trefftz FESpace.
  template <typename SCAL>
  void
  etFesVTransformV (SliceVector<SCAL> &vec, const TRANSFORM_TYPE type,
                    const Matrix<SCAL> &elmat, const Matrix<SCAL> &elmat_inv)
  {
    const size_t ndof = elmat.Width ();

    if (type == TRANSFORM_RHS)
      {
        Vector<SCAL> new_vec (ndof);
        new_vec = Trans (elmat) * vec;
        vec.Range (ndof) = new_vec;
      }
    else if (type == TRANSFORM_SOL)
      {
        Vector<SCAL> new_vec (vec.Size ());
        new_vec = 0;
        new_vec.Range (0, elmat.Height ()) = elmat * vec.Range (0, ndof);
        vec = new_vec;
      }
    else if (type == TRANSFORM_SOL_INVERSE)
      {
        if (vec.Size () != elmat.Height ())
          throw std::invalid_argument (
              "given vec does not match the needed dimension.");

        // const auto elmat_inv = getPseudoInverse (elmat, 0);
        Vector<SCAL> tmp_vec (elmat_inv.Height ());
        tmp_vec = elmat_inv * vec;
        vec = tmp_vec;
      }
    else
      {
        stringstream err;
        err << "VTransformV is not implemented for TRANSFORM_TYPE " << type;
        throw std::invalid_argument (err.str ());
      }
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformVR (ElementId ei, SliceVector<double> vec,
                                      TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: VTransform");
    RegionTimer reg (timer);

    etFesVTransformV (vec, type, *this->etmats[ei.Nr ()],
                      *this->etmats_inv[ei.Nr ()]);
  }

  template <typename T>
  void
  EmbTrefftzFESpace<T>::VTransformVC (ElementId ei, SliceVector<Complex> vec,
                                      TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: VTransform");
    RegionTimer reg (timer);

    etFesVTransformV (vec, type, *this->etmatsc[ei.Nr ()],
                      *this->etmatsc_inv[ei.Nr ()]);
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

template <typename T> string EmbTrefftzFESpace<T>::GetClassName () const
{
  return this->name;
}

////////////////////////// python interface ///////////////////////////

#ifdef NGS_PYTHON
template <typename T> void ExportETSpace (py::module m, string label)
{
  auto pyspace
      = ngcomp::ExportFESpace<ngcomp::EmbTrefftzFESpace<T>> (m, label);

  // pyspace.def (py::init ([pyspace] (shared_ptr<T> fes) {
  // py::list info;
  // auto ma = fes->GetMeshAccess ();
  // info.append (ma);
  // auto nfes = make_shared<ngcomp::EmbTrefftzFESpace<T>> (fes);
  // nfes->Update ();
  // nfes->FinalizeUpdate ();
  // connect_auto_update (nfes.get ());
  // return nfes;
  //}),
  // py::arg ("fes"));
}

void ExportEmbTrefftz (py::module m)
{
  ExportETSpace<ngcomp::L2HighOrderFESpace> (m, "L2EmbTrefftzFESpace");
  ExportETSpace<ngcomp::VectorL2FESpace> (m, "VectorL2EmbTrefftzFESpace");
  ExportETSpace<ngcomp::MonomialFESpace> (m, "MonomialEmbTrefftzFESpace");
  ExportETSpace<ngcomp::CompoundFESpace> (m, "CompoundEmbTrefftzFESpace");

  m.def (
      "EmbeddedTrefftzFES",
      [] (shared_ptr<ngcomp::TrefftzEmbedding> emb)
          -> shared_ptr<ngcomp::FESpace> {
        shared_ptr<ngcomp::FESpace> fes = emb->GetFES ();
        shared_ptr<ngcomp::FESpace> nfes;
        if (dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::L2HighOrderFESpace>> (emb);
        else if (dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::VectorL2FESpace>> (emb);
        else if (dynamic_pointer_cast<ngcomp::MonomialFESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::MonomialFESpace>> (emb);
        else if (dynamic_pointer_cast<ngcomp::CompoundFESpace> (fes))
          nfes = make_shared<
              ngcomp::EmbTrefftzFESpace<ngcomp::CompoundFESpace>> (emb);
        else
          throw Exception ("Unknown base fes");
        return nfes;
      },
      R"mydelimiter(
        Given a TrefftzEmbedding this wrapper produces a Trefftz FESpace using local projections,
        following the Embedded Trefftz-DG methodology.

        :param TrefftzEmbedding: The Trefftz embedding object.

        :return: EmbTrefftzFES
        )mydelimiter",
      py::arg ("emb"));

  py::class_<TrefftzEmbedding, shared_ptr<TrefftzEmbedding>> (
      m, "TrefftzEmbedding",
      R"mydelimiter(
                Gives access to the embedding matrix and a particular solution and
                can be used to construct an EmbeddedTrefftzFESpace.

                The dimension of the local Trefftz space is determined by the kernel of `top`,
                after removing the dofs fixed by the conforming condition in `cop` and `crhs`.

                If a different test space is used, the dimension of the local Trefftz space is
                at best dim(fes)-dim(fes_test) and may increase by zero singular values of `top`
                (with respect to the threshold `eps`).
            )mydelimiter")
      .def (py::init<shared_ptr<SumOfIntegrals>, shared_ptr<SumOfIntegrals>,
                     shared_ptr<SumOfIntegrals>, shared_ptr<SumOfIntegrals>,
                     size_t, double, shared_ptr<FESpace>, shared_ptr<FESpace>,
                     shared_ptr<FESpace>, shared_ptr<BitArray>> (),
            R"mydelimiter(
                Constructs a new Trefftz embedding object.

                 :param top: the differential operation. Can be None
                 :param trhs: right hand side of the var. formulation
                 :param cop: left hand side of the conformity operation
                 :param crhs: right hand side of the conformity operation
                 :param eps: cutoff for singular values from the SVD of the local operator.
                        values below eps are considered zero and therefore in the kernel of `top`.
                        (default: 0.0)
                 :param ndof_trefftz: fixes the number of degrees of freedom per element
                     that are to be considered in the Trefftz space generated by `top`
                     (i.e. the local dimension of the kernel of `top` on one element)
                     cannot be used together with `eps` (default: 0)
                 :param fes: the finite element space of `top` (optional, determined
                     from `top` if not given)
                 :param fes_test: the finite element test space of `top` (optional,
                     determined from `top` if not given)
                 :param fes_conformity: finite element space of the conformity operation (optional,
                     determined from `cop` if not given)
            )mydelimiter",
            py::arg ("top") = nullptr, py::arg ("trhs") = nullptr,
            py::arg ("cop") = nullptr, py::arg ("crhs") = nullptr,
            py::arg ("ndof_trefftz") = 0, py::arg ("eps") = 0.0,
            py::arg ("fes") = nullptr, py::arg ("fes_test") = nullptr,
            py::arg ("fes_conformity") = nullptr,
            py::arg ("ignoredofs") = nullptr)
      .def ("Embed", &ngcomp::TrefftzEmbedding::Embed,
            "Embed a Trefftz GridFunction into the underlying FESpace")
      .def ("GetEmbedding", &ngcomp::TrefftzEmbedding::GetEmbedding,
            "Get the sparse embedding matrix")
      .def ("GetParticularSolution",
            static_cast<shared_ptr<const BaseVector> (TrefftzEmbedding::*) ()
                            const> (
                &ngcomp::TrefftzEmbedding::GetParticularSolution),
            "Particular solution as GridFunction vector of the underlying "
            "FESpace")
      .def ("GetParticularSolution",
            static_cast<shared_ptr<const BaseVector> (TrefftzEmbedding::*) (
                shared_ptr<SumOfIntegrals>) const> (
                &ngcomp::TrefftzEmbedding::GetParticularSolution),
            "Particular solution as GridFunction vector of the underlying "
            "FESpace, given a trhs");
}

#endif // NGS_PYTHON
