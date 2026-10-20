#include "embtrefftz.hpp"
#include "monomialfespace.hpp"
#include "trefftzfespace.hpp"
#include "lapack_loader.hpp"
#ifdef NGSTREFFTZ_USE_XFEM
#include <xfem/aggregates.hpp>
#endif
#include <cfloat>
#include <compressedfespace.hpp>
#include <type_traits>

using namespace ngbla;
using namespace ngcomp;

namespace
{
  shared_ptr<DofPatches> GetDofPatches (const FESpace &space)
  {
    if (auto provider = dynamic_cast<const DofPatchProvider *> (&space))
      return make_shared<DofPatches> (provider->GetDofPatches ());
    if (auto compressed = dynamic_cast<const CompressedFESpace *> (&space))
      {
        auto patches = GetDofPatches (*compressed->GetBaseSpace ());
        if (patches)
          {
            auto flat_dofs = patches->dofs.AsArray ();
            Array<DofId> dofs (flat_dofs.Size (), flat_dofs.Data ());
            compressed->WrapDofs (dofs);
          }
        return patches;
      }
    return nullptr;
  }

  /// @returns `max(a-b, 0)` without integer underflow risk, if `b > a`.
  INLINE size_t max_a_minus_b_or_0 (const size_t a, const size_t b) noexcept
  {
    return (a >= b) ? a - b : 0;
  }

  template <typename SCAL, typename TDIST>
  inline void AddIntegrationToElementMatrix (
      MatrixView<SCAL, RowMajor, size_t, size_t, TDIST> elmat,
      const Array<shared_ptr<BilinearFormIntegrator>> &bf_integrators,
      const MeshAccess &ma, const ElementId &element_id, const FESpace &fes,
      const FESpace &fes_test, LocalHeap &lh)
  {
    const HeapReset hr (lh);
    if (!elmat.Height () || !elmat.Width () || !bf_integrators.Size ())
      return;

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
        if (bfi->DefinedOnElement (element_id.Nr ())
            && (!bfi->GetDefinedOn ().Size ()
                || bfi->DefinedOn (ma.GetElIndex (element_id))))
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

  INLINE bool IsInactiveOrHiddenDof (const FESpace &fes, DofId d)
  {
    if (!IsRegularDof (d))
      return true;
    return fes.GetDofCouplingType (d) <= HIDDEN_DOF;
  }

  INLINE bool
  IsIgnoredRegularDof (shared_ptr<const BitArray> ignoredofs, DofId d)
  {
    if (!ignoredofs)
      return false;
    if (!IsRegularDof (d))
      return false;
    const size_t id = size_t (d);
    return (id < ignoredofs->Size ()) && ignoredofs->Test (id);
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
      const FESpace &fes, const FESpace &fes_test,
      const shared_ptr<const FESpace> fes_conformity, Array<DofId> &dofs,
      Array<DofId> &dofs_test, Array<DofId> &conformity_dofs, LocalHeap &lh,
      bool compute_new_dofs = false, shared_ptr<BitArray> ignoredofs = nullptr)
  {
    static Timer timer ("EmbTrefftz: extractVisibleDofs");
    RegionTimer reg (timer);
    Array<int> pos_trial, pos_test, pos_conf;
    Array<DofId> vdofs, vdofs_test, vdofs_conformity;

    pos_trial.SetSize0 ();
    pos_test.SetSize0 ();
    pos_conf.SetSize0 ();

    vdofs.SetSize0 ();
    vdofs_test.SetSize0 ();
    vdofs_conformity.SetSize0 ();

    for (size_t j = 0; j < dofs.Size (); j++)
      {
        DofId d = dofs[j];
        if (IsInactiveOrHiddenDof (fes, d))
          continue;
        if (IsIgnoredRegularDof (ignoredofs, d))
          continue;
        pos_trial.Append (j);
        vdofs.Append (d);
      }

    for (size_t i = 0; i < dofs_test.Size (); i++)
      {
        DofId d = dofs_test[i];
        if (IsInactiveOrHiddenDof (fes_test, d))
          continue;
        pos_test.Append (i);
        vdofs_test.Append (d);
      }

    if (fes_conformity)
      for (size_t i = 0; i < conformity_dofs.Size (); i++)
        {
          DofId d = conformity_dofs[i];
          if (IsInactiveOrHiddenDof (*fes_conformity, d))
            continue;
          pos_conf.Append (i);
          vdofs_conformity.Append (d);
        }

    FlatMatrix<SCAL> velmat (pos_test.Size () + pos_conf.Size (),
                             pos_trial.Size (), lh);

    size_t conformity_offset_old
        = fes_conformity ? conformity_dofs.Size () : 0;
    size_t conformity_offset_new = fes_conformity ? pos_conf.Size () : 0;

    if (fes_conformity)
      for (size_t vi = 0; vi < pos_conf.Size (); vi++)
        for (size_t vj = 0; vj < pos_trial.Size (); vj++)
          velmat (vi, vj) = elmat (pos_conf[vi], pos_trial[vj]);

    for (size_t vi = 0; vi < pos_test.Size (); vi++)
      for (size_t vj = 0; vj < pos_trial.Size (); vj++)
        velmat (conformity_offset_new + vi, vj)
            = elmat (conformity_offset_old + pos_test[vi], pos_trial[vj]);

    if (compute_new_dofs)
      {
        conformity_dofs = std::move (vdofs_conformity);
        dofs = std::move (vdofs);
        dofs_test = std::move (vdofs_test);
      }

    return velmat;
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
          throw Exception (
              "Two different trial FESpaces in the same operator");
        else if (proxy && proxy->IsTestFunction () && !fes_test)
          fes_test = proxy->GetFESpace ();
        else if (proxy && proxy->IsTestFunction ()
                 && fes_test != proxy->GetFESpace ())
          throw Exception ("Two different test FESpaces in the same operator");
      });
  }

  // Global numbering is [conformity DOFs | local Trefftz DOFs | ignored DOFs].
  Array<DofId>
  NumberIgnoredDofs (const FESpace &fes, shared_ptr<const BitArray> ignoredofs,
                     size_t offset)
  {
    Array<DofId> numbers;
    if (ignoredofs)
      {
        numbers.SetSize (fes.GetNDof ());
        numbers = NO_DOF_NR;
        for (size_t d : Range (ignoredofs->Size ()))
          if (ignoredofs->Test (d))
            numbers[d] = offset++;
      }
    return numbers;
  }

  /// Number element columns without assembling a global embedding matrix.
  template <typename SCAL>
  Table<DofId>
  ElementEmbeddingDofs (const Array<optional<Matrix<SCAL>>> &etmats,
                        FlatArray<size_t> local_ndofs_trefftz,
                        const FESpace &fes,
                        shared_ptr<const FESpace> fes_conformity,
                        shared_ptr<const BitArray> ignoredofs)
  {
    auto ma = fes.GetMeshAccess ();
    VorB vb = ma->GetNE (VOL) ? VOL : BND;
    size_t nc = fes_conformity ? fes_conformity->GetNDof () : 0;
    size_t nt = 0;
    for (size_t n : local_ndofs_trefftz)
      nt += n;
    auto ignored_numbers = NumberIgnoredDofs (fes, ignoredofs, nc + nt);
    TableCreator<DofId> cols (ma->GetNE (vb));
    for (; !cols.Done (); cols++)
      {
        size_t next = nc;
        for (auto ei : ma->Elements (vb))
          if (etmats[ei.Nr ()])
            {
              Array<DofId> dofs;
              if (fes_conformity)
                {
                  fes_conformity->GetDofNrs (ei, dofs, VISIBLE_DOF);
                  for (DofId d : dofs)
                    if (IsRegularDof (d))
                      cols.Add (ei.Nr (), d);
                }
              for (size_t j : Range (local_ndofs_trefftz[ei.Nr ()]))
                cols.Add (ei.Nr (), next + j);
              next += local_ndofs_trefftz[ei.Nr ()];
              if (ignoredofs)
                {
                  fes.GetDofNrs (ei, dofs);
                  for (DofId d : dofs)
                    if (IsIgnoredRegularDof (ignoredofs, d))
                      cols.Add (ei.Nr (), ignored_numbers[d]);
                }
            }
      }
    return cols.MoveTable ();
  }

  template <typename SCAL>
  void
  AssemblePatchEmbeddings (const Array<optional<Matrix<SCAL>>> &local_matrices,
                           FlatArray<LocalEmbeddingData> local_data,
                           const DofPatches &patches, const FESpace &fes,
                           shared_ptr<const FESpace> fes_conformity,
                           shared_ptr<const BitArray> ignoredofs,
                           Array<optional<Matrix<SCAL>>> &etmats,
                           Table<DofId> &tdof_nrs,
                           Array<size_t> &local_ndofs_trefftz)
  {
    size_t nc = fes_conformity ? fes_conformity->GetNDof () : 0;
    size_t nt = 0;
    for (const auto &data : local_data)
      nt += data.ndof_trefftz;
    auto ignored_numbers = NumberIgnoredDofs (fes, ignoredofs, nc + nt);

    auto ma = fes.GetMeshAccess ();
    VorB vb = ma->GetNE (VOL) ? VOL : BND;
    Array<int> element_patch (ma->GetNE (vb));
    element_patch = -1;
    for (size_t p : Range (local_data))
      if (local_data[p].active && local_data[p].trial_dofs.Size ())
        for (size_t el : patches.elements[p])
          element_patch[el] = p;

    TableCreator<DofId> cols (ma->GetNE (vb));
    for (; !cols.Done (); cols++)
      for (auto ei : ma->Elements (vb))
        {
          Array<DofId> dofs;
          fes.GetDofNrs (ei, dofs);
          if (!dofs.Size ())
            continue;
          int p = element_patch[ei.Nr ()];
          if (p >= 0)
            {
              const auto &data = local_data[p];
              for (DofId d : data.conformity_dofs)
                cols.Add (ei.Nr (), d);
              for (size_t j : Range (data.ndof_trefftz))
                cols.Add (ei.Nr (), nc + data.trefftz_dof_offset + j);
            }
          for (DofId d : dofs)
            if (IsIgnoredRegularDof (ignoredofs, d))
              cols.Add (ei.Nr (), ignored_numbers[d]);
        }
    tdof_nrs = cols.MoveTable ();
    etmats.SetSize (ma->GetNE (vb));
    local_ndofs_trefftz.SetSize (ma->GetNE (vb));
    local_ndofs_trefftz = 0;
    for (auto ei : ma->Elements (vb))
      if (tdof_nrs[ei.Nr ()].Size ())
        {
          Array<DofId> dofs;
          fes.GetDofNrs (ei, dofs);
          auto tdofs = tdof_nrs[ei.Nr ()];
          auto &elmat = etmats[ei.Nr ()].emplace (dofs.Size (), tdofs.Size ());
          elmat = SCAL (0.0);
          int p = element_patch[ei.Nr ()];
          if (p >= 0)
            {
              const auto &data = local_data[p];
              const auto &matrix = *local_matrices[p];
              for (size_t i : Range (dofs))
                {
                  size_t row = data.trial_dofs.Pos (dofs[i]);
                  if (row != Array<DofId>::ILLEGAL_POSITION)
                    elmat.Row (i).Range (matrix.Width ()) = matrix.Row (row);
                }
              local_ndofs_trefftz[ei.Nr ()] = data.ndof_trefftz;
            }
          for (size_t i : Range (dofs))
            if (IsIgnoredRegularDof (ignoredofs, dofs[i]))
              elmat (i, tdofs.Pos (ignored_numbers[dofs[i]])) = SCAL (1.0);
        }
  }

  /// The sparse matrix is an optional view of the stored element matrices.
  template <typename SCAL>
  shared_ptr<BaseMatrix>
  Elmats2Sparse (const Array<optional<Matrix<SCAL>>> &etmats,
                 const Table<DofId> &tdof_nrs, const FESpace &fes, size_t ndof,
                 shared_ptr<const BitArray> ignoredofs, bool patchwise)
  {
    auto ma = fes.GetMeshAccess ();
    VorB vb = ma->GetNE (VOL) ? VOL : BND;
    TableCreator<DofId> rows (ma->GetNE (vb));
    for (; !rows.Done (); rows++)
      for (auto ei : ma->Elements (vb))
        if (etmats[ei.Nr ()])
          {
            Array<DofId> dofs;
            fes.GetDofNrs (ei, dofs);
            for (DofId d : dofs)
              if (IsRegularDof (d))
                rows.Add (ei.Nr (), d);
          }
    auto row_table = rows.MoveTable ();
    Table<DofId> graph (tdof_nrs);
    auto matrix = make_shared<SparseMatrix<SCAL>> (fes.GetNDof (), ndof,
                                                   row_table, graph, false);
    matrix->SetZero ();
    for (auto ei : ma->Elements (vb))
      if (etmats[ei.Nr ()])
        {
          Array<DofId> dofs;
          fes.GetDofNrs (ei, dofs);
          auto cols = tdof_nrs[ei.Nr ()];
          const auto &elmat = *etmats[ei.Nr ()];
          for (size_t i : Range (dofs))
            if (IsRegularDof (dofs[i]))
              {
                for (size_t j : Range (cols))
                  if (patchwise)
                    (*matrix) (dofs[i], cols[j]) = elmat (i, j);
                  else
                    (*matrix) (dofs[i], cols[j]) += elmat (i, j);
              }
        }
    if (ignoredofs)
      {
        size_t next = ndof - ignoredofs->NumSet ();
        for (size_t d : Range (ignoredofs->Size ()))
          if (ignoredofs->Test (d))
            (*matrix) (d, next++) = SCAL (1.0);
      }
    return matrix;
  }

  void calculateBilinearFormIntegrators (
      const SumOfIntegrals &bf,
      Array<shared_ptr<BilinearFormIntegrator>> bfis[4],
      shared_ptr<MeshAccess> ma)
  {
    for (auto icf : bf.icfs)
      {
        DifferentialSymbol &dx = icf->dx;
        shared_ptr<BilinearFormIntegrator> bfi
            = icf->MakeBilinearFormIntegrator ();

        if (dx.definedon) // DifferentialSymbol cannot set definedon region as
                          // it has no mesh access.
          {
            if (auto definedon_string = get_if<string> (&*dx.definedon);
                definedon_string)
              {
                Region reg (ma, dx.vb, *definedon_string);
                bfi->SetDefinedOn (reg.Mask ());
              }
          }

        bfis[dx.vb] += bfi;
      }
  }

  void CalculateLinearFormIntegrators (
      const SumOfIntegrals &trhs,
      Array<shared_ptr<LinearFormIntegrator>> lfis[4],
      shared_ptr<MeshAccess> ma, bool volume_only = false)
  {
    for (auto icf : trhs.icfs)
      {
        DifferentialSymbol &dx = icf->dx;
        if (volume_only && dx.vb != VOL)
          throw Exception ("Patch right-hand sides require volume integrals");
        shared_ptr<LinearFormIntegrator> lfi
            = icf->MakeLinearFormIntegrator ();

        if (dx.definedon)
          {
            if (auto definedon_string = get_if<string> (&*dx.definedon);
                definedon_string)
              {
                Region reg (ma, dx.vb, *definedon_string);
                lfi->SetDefinedOn (reg.Mask ());
              }
          }

        lfis[dx.vb] += lfi;
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
                             const Ngs_Element &mesh_element,
                             shared_ptr<MeshAccess> ma)
  {
    for (auto icf : bf.icfs)
      {
        if (icf->dx.vb != mesh_element.VB ())
          continue;

        if (icf->dx.definedonelements
            && !icf->dx.definedonelements->Test (mesh_element.Nr ()))
          continue;

        if (icf->dx.definedon) // check region
          {
            if (auto definedon_string = get_if<string> (&*icf->dx.definedon);
                definedon_string)
              {
                Region reg (ma, icf->dx.vb, *definedon_string);
                const BitArray &definedon = reg.Mask ();
                size_t mat = ma->GetElIndex (mesh_element);
                if (mat >= definedon.Size ())
                  continue;
                if (definedon.Test (mat))
                  return true;
                continue;
              }
          }

        return true;
      }
    return false;
  }
}

namespace
{
  void
  LapackSVD (SliceMatrix<double, ColMajor> A, SliceMatrix<double, ColMajor> U,
             SliceMatrix<double, ColMajor> V)
  {
    static Timer t ("EmbTrefftz: LapackSVD");
    RegionTimer reg (t);
    ngbla::integer n = A.Width (), m = A.Height ();
    Vector<> S (min (n, m));
    Array<double> work (4 * m * m + 6 * m + m + 100);
    Array<ngbla::integer> iwork (max (n, m) * 9);
    ngbla::integer info;
    char jobz = 'A';
    ngbla::integer lda = A.Dist (), ldu = U.Dist (), ldv = V.Dist ();
    ngbla::integer lwork = work.Size ();

    auto dgesdd = ngstrefftz::GetLapackGesdd ().realgesdd;
    if (!dgesdd)
      throw Exception ("LAPACK dgesdd symbol is not available");
    dgesdd (&jobz, &m, &n, A.Data (), &lda, S.Data (), U.Data (), &ldu,
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
    static Timer t ("EmbTrefftz: LapackSVD");
    RegionTimer reg (t);
    ngbla::integer n = A.Width (), m = A.Height ();
    Vector<> S (min (n, m));
    Array<Complex> work (4 * m * m + 6 * m + m + 100);
    Array<ngbla::integer> iwork (max (n, m) * 9);
    Array<double> rwork (5 * max (n, m) * min (n, m) + 5 * min (n, m));
    ngbla::integer info;
    char jobz = 'A';
    ngbla::integer lda = A.Dist (), ldu = U.Dist (), ldv = V.Dist ();
    ngbla::integer lwork = work.Size ();

    auto zgesdd = ngstrefftz::GetLapackGesdd ().complexgesdd;
    if (!zgesdd)
      throw Exception ("LAPACK zgesdd symbol is not available");
    zgesdd (&jobz, &m, &n, A.Data (), &lda, S.Data (), U.Data (), &ldu,
            V.Data (), &ldv, work.Data (), &lwork, rwork.Data (),
            iwork.Data (), &info);
    if (info != 0)
      throw Exception ("something went wrong in the svd "
                       + std::to_string (info));
    A = 0.0;
    A.Diag (0) = S;
  }

  /// Returns V and UT = U^H from A = U * Sigma * V^H.
  /// A gets overwritten with Sigma.
  /// @param A has dimension n * m
  /// @param U has dimension n * n
  /// @param V has dimension m * m
  template <typename SCAL, typename TDIST>
  void
  GetSVD (MatrixView<SCAL, ngbla::RowMajor, size_t, size_t, TDIST> A_to_sigmaT,
          FlatMatrix<SCAL, ColMajor> UT, FlatMatrix<SCAL, ColMajor> V)
  {
    auto [Aheight, Awidth] = A_to_sigmaT.Shape ();
    FlatMatrix<SCAL, ColMajor> AT (Awidth, Aheight, &A_to_sigmaT (0, 0));
    // FlatMatrix<SCAL, ColMajor> VT (Aheight, Aheight, &UT(0,0));
    // FlatMatrix<SCAL, ColMajor> U (Awidth, Awidth, &V(0,0));
    // Matrix<SCAL, ColMajor> AA (height, width);
    // AA = A;

    const auto &lapack = ngstrefftz::GetLapackGesdd ();
    if constexpr (std::is_same_v<SCAL, double>)
      {
        if (lapack.realgesdd)
          LapackSVD (AT, V, UT);
        else
          CalcSVD (AT, V, UT);
      }
    else if constexpr (std::is_same_v<SCAL, Complex>)
      {
        if (lapack.complexgesdd)
          LapackSVD (AT, V, UT);
        else
          throw Exception ("complex SVD requires LAPACK zgesdd");
        // LAPACK factors A^T; conjugate both factors to obtain those of A.
        V = Conj (V);
        UT = Conj (UT);
      }
    else
      throw Exception ("SVD only implemented for double and Complex");

    // A = static_cast<SCAL> (0.0);
    //  A.Diag(0)=AA.Diag();
    // for (size_t i = 0; i < min (A.Width (), A.Height ()); i++)
    // A (i, i) = AA (i, i);
  }
}

namespace
{
  /// @param num_zeros_guess Guess (lower bound) for the number of singular
  ///     values `sigma_i = 0` in the matrix Sigma. Is useful when there are
  ///     very-close-to-zero singular values that you can detect in advance.
  /// @return pseudoinverse of A (as some `ngbla::Expr` type to avoid
  /// allocations)
  template <typename SCAL, typename TDIST, ORDERING SIG_ORD>
  Matrix<SCAL>
  invertSVD (const FlatMatrix<SCAL, ColMajor> &UT,
             const MatrixView<SCAL, SIG_ORD, size_t, size_t, TDIST> &Sigma,
             const FlatMatrix<SCAL, ColMajor> &V, size_t num_zeros_guess,
             LocalHeap &lh, bool take_sqrt = false)
  {
    const HeapReset hr (lh);
    const auto [m, n] = Sigma.Shape ();
    const size_t diag_len = min (m, n);

    // find out if the guess was too low
    size_t num_zeros = num_zeros_guess;
    for (size_t i = num_zeros_guess + 1;
         i <= diag_len && Sigma (diag_len - i, diag_len - i) == 0.; i++)
      num_zeros = i;

    FlatMatrix<SCAL> sigma_inv_times_ut (n, UT.Width (), lh);
    const size_t nonzero_diag_len = min (n - num_zeros, m);
    for (size_t i = 0; i < nonzero_diag_len; i++)
      if (take_sqrt)
        sigma_inv_times_ut.Row (i) = 1.0 / sqrt (Sigma (i, i)) * UT.Row (i);
      else
        sigma_inv_times_ut.Row (i) = 1.0 / Sigma (i, i) * UT.Row (i);
    sigma_inv_times_ut.Rows (nonzero_diag_len, n) = SCAL (0.);

    return V * sigma_inv_times_ut;
  }

  /// @returns the pseudo-inverse of `mat`
  /// @param num_zero number of (near) zero singular values in `mat`
  template <typename SCAL>
  Matrix<SCAL>
  getPseudoInverse (const FlatMatrix<SCAL> mat, const size_t num_zero)
  {
    LocalHeap lh = LocalHeap (100 * 1000 * 1000, "embt");
    const auto [n, m] = mat.Shape ();
    FlatMatrix<SCAL> sigma (n, m, lh);
    FlatMatrix<SCAL, ColMajor> UT (n, n, lh);
    FlatMatrix<SCAL, ColMajor> V (m, m, lh);

    Matrix<SCAL> elmat_inv (m, n);

    sigma = mat;
    GetSVD (sigma, UT, V);
    elmat_inv = invertSVD (UT, sigma, V, num_zero, lh);
    return elmat_inv;
  }

  template <typename SCAL>
  void getPseudoInverse (FlatMatrix<SCAL> mat, const size_t num_zero,
                         LocalHeap &lh, int take_sqrt = false)
  {
    HeapReset hr (lh);
    const auto [n, m] = mat.Shape ();
    FlatMatrix<SCAL> sigma (n, m, lh);
    FlatMatrix<SCAL, ColMajor> UT (n, n, lh);
    FlatMatrix<SCAL, ColMajor> V (m, m, lh);

    sigma = mat;
    GetSVD (sigma, UT, V);
    mat = invertSVD (UT, sigma, V, num_zero, lh, take_sqrt);
  }

  template <typename SCAL>
  FlatVector<SCAL>
  AssembleElementRhs (Array<shared_ptr<LinearFormIntegrator>> lfis[4],
                      const FESpace &fes_test, const ElementId ei,
                      const MeshAccess &ma, bool has_inactive_dofs,
                      LocalHeap &mlh)
  {
    Array<DofId> dofs_test;
    fes_test.GetDofNrs (ei, dofs_test);

    const auto &test_fel = fes_test.GetFE (ei, mlh);
    const auto &trafo = ma.GetTrafo (ei, mlh);
    FlatVector<SCAL> elvec (dofs_test.Size (), mlh),
        elveci (dofs_test.Size (), mlh);
    elvec = static_cast<SCAL> (0.0);

    for (const auto vorb : { VOL, BND, BBND, BBBND })
      {
        for (const auto &lfi : lfis[vorb])
          {
            if (lfi->DefinedOnElement (ei.Nr ())
                && lfi->DefinedOn (ma.GetElIndex (ei)))
              {
                auto &mapped_trafo = trafo.AddDeformation (
                    lfi->GetDeformation ().get (), mlh);
                lfi->CalcElementVector (test_fel, mapped_trafo, elveci, mlh);
                elvec += elveci;
              }
          }
      }

    if (has_inactive_dofs)
      {
        Array<int> pos_test;
        pos_test.SetSize0 ();
        for (size_t i = 0; i < dofs_test.Size (); i++)
          if (!IsInactiveOrHiddenDof (fes_test, dofs_test[i]))
            pos_test.Append (i);

        FlatVector<SCAL> velvec (pos_test.Size (), mlh);
        for (size_t vi = 0; vi < pos_test.Size (); vi++)
          velvec[vi] = elvec[pos_test[vi]];

        elvec.AssignMemory (velvec.Size (), &velvec[0]);
      }
    return elvec;
  }

}

namespace ngcomp
{
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

  template size_t
  ngcomp::calcNdofTrefftz<double> (size_t, size_t, size_t,
                                   std::variant<size_t, double>, bool,
                                   ngbla::SliceVector<double>);

  template size_t
  ngcomp::calcNdofTrefftz<Complex> (size_t, size_t, size_t,
                                    std::variant<size_t, double>, bool,
                                    ngbla::SliceVector<Complex>);

  template <typename SCAL> struct TrefftzEmbedding::EmbeddingWork
  {
    Array<shared_ptr<BilinearFormIntegrator>> op_integrators[4],
        cop_lhs_integrators[4], cop_rhs_integrators[4], fes_ip_integrators[4];
    Array<shared_ptr<LinearFormIntegrator>> lfis[4];
    Array<optional<Matrix<SCAL>>> local_matrices, inverses;
    shared_ptr<VVector<SCAL>> particular_solution;
    bool any_fes_has_inactive_dofs = false;
    bool test_has_inactive_dofs = false;
    bool patchwise = false;

    // statistics stuff
    Vector<double> sing_val_avg, sing_val_max, sing_val_min, sing_val_count;
    double conformity_residual = 0.0, trefftz_residual = 0.0;
    mutex stats_mutex;
  };

  template <typename SCAL, typename GETDOFS, typename ADDMATRIX,
            typename ADDVECTOR>
  INLINE void TrefftzEmbedding::ComputeLocalEmbedding (
      size_t nr, EmbeddingWork<SCAL> &work, const GETDOFS &get_dofs,
      const ADDMATRIX &add_matrix, const ADDVECTOR &add_vector, LocalHeap &lh)
  {
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

    LocalEmbeddingData data;
    auto &dofs = data.trial_dofs;
    auto &dofs_test = data.test_dofs;
    auto &dofs_conforming = data.conformity_dofs;
    if (work.patchwise)
      dofs = patches->dofs[nr];
    else
      get_dofs (*fes, dofs);
    if (fes_test == fes)
      dofs_test = dofs;
    else
      get_dofs (*fes_test, dofs_test);
    if (fes_conformity)
      get_dofs (*fes_conformity, dofs_conforming);
    if (!dofs.Size ())
      return;

    size_t ndof = dofs.Size ();
    size_t ndof_test = dofs_test.Size ();
    size_t ndof_conforming = dofs_conforming.Size ();

    FlatMatrix<SCAL> elmat_A (ndof_test + ndof_conforming, ndof, lh);
    FlatMatrix<SCAL> elmat_B (ndof_test + ndof_conforming, ndof_conforming,
                              lh);
    Matrix<SCAL> elmat_A_copy;
    elmat_A = static_cast<SCAL> (0.);
    elmat_B = static_cast<SCAL> (0.);

    // split conformity part
    auto [elmat_Cl, elmat_L] = elmat_A.SplitRows (ndof_conforming);
    auto elmat_Cr = elmat_B.Rows (ndof_conforming);

    // the diff. operator L operates only on volume/bnd (vb) terms
    add_matrix (elmat_L, work.op_integrators[vb], *fes, *fes_test, dofs,
                dofs_test, lh);
    if (fes_conformity)
      {
        for (const auto vorb : { VOL, BND, BBND, BBBND })
          {
            add_matrix (elmat_Cl, work.cop_lhs_integrators[vorb], *fes,
                        *fes_conformity, dofs, dofs_conforming, lh);
            add_matrix (elmat_Cr, work.cop_rhs_integrators[vorb],
                        *fes_conformity, *fes_conformity, dofs_conforming,
                        dofs_conforming, lh);
          }
      }

    if (stats)
      elmat_A_copy = elmat_A;

    FlatMatrix<SCAL> fes_ip_sqinv;
    Array<DofId> weight_dofs;
    if (fes_ip)
      {
        if (work.any_fes_has_inactive_dofs || ignoredofs)
          weight_dofs = dofs;
        fes_ip_sqinv.AssignMemory (ndof, ndof, lh);
        fes_ip_sqinv = 0.;
        add_matrix (fes_ip_sqinv, work.fes_ip_integrators[vb], *fes, *fes,
                    dofs, dofs, lh);
      }

    if (work.any_fes_has_inactive_dofs || ignoredofs)
      {
        if (stats)
          elmat_A_copy = extractVisibleDofs (
              elmat_A_copy, *fes, *fes_test, fes_conformity, dofs, dofs_test,
              dofs_conforming, lh, false, ignoredofs);
        if (fes_conformity)
          {
            elmat_B.Assign (extractVisibleDofs (
                elmat_B, *fes_conformity, *fes_test, fes_conformity,
                dofs_conforming, dofs_test, dofs_conforming, lh, false));
          }

        elmat_A.Assign (extractVisibleDofs (
            elmat_A, *fes, *fes_test, fes_conformity, dofs, dofs_test,
            dofs_conforming, lh, true, ignoredofs));

        ndof = dofs.Size ();
        ndof_test = dofs_test.Size ();
        ndof_conforming = dofs_conforming.Size ();
      }
    if (fes_ip && ndof)
      {
        if (ndof != fes_ip_sqinv.Height ())
          {
            FlatMatrix<SCAL> weight (ndof, ndof, lh);
            Array<size_t> positions (ndof);
            for (size_t i : Range (dofs))
              positions[i] = weight_dofs.Pos (dofs[i]);
            for (size_t i : Range (dofs))
              for (size_t j : Range (dofs))
                weight (i, j) = fes_ip_sqinv (positions[i], positions[j]);
            fes_ip_sqinv.Assign (weight);
          }
        getPseudoInverse (fes_ip_sqinv, 0, lh, true);
        FlatMatrix<SCAL> weighted_A (elmat_A.Height (), ndof, lh);
        weighted_A = elmat_A * fes_ip_sqinv;
        elmat_A.Assign (weighted_A);
      }

    Matrix<SCAL> elmat_T;
    size_t ndof_trefftz_i = ndof;
    if (!ndof || !(ndof_test + ndof_conforming))
      {
        elmat_T.SetSize (ndof, ndof_conforming + ndof);
        elmat_T = SCAL (0.0);
        if (ndof)
          {
            if (fes_ip)
              elmat_T = fes_ip_sqinv;
            else
              elmat_T = Identity (ndof);
          }
        work.inverses[nr].emplace (ndof, ndof_test);
        *work.inverses[nr] = SCAL (0.0);
      }
    else
      {
        FlatMatrix<SCAL, ColMajor> UT (elmat_A.Height (), lh),
            V (elmat_A.Width (), lh);
        GetSVD<SCAL> (elmat_A, UT, V);

        ndof_trefftz_i
            = calcNdofTrefftz (ndof, ndof_test, ndof_conforming, ndof_trefftz,
                               !top, elmat_A.Diag (0));
        if (ndof_trefftz_i + ndof_conforming == 0)
          throw std::invalid_argument ("zero trefftz dofs");

        Matrix<SCAL> elmat_A_inv
            = invertSVD (UT, elmat_A, V, ndof_trefftz_i, lh);

        // T = (T_c | T_t)
        elmat_T.SetSize (ndof, ndof_conforming + ndof_trefftz_i);
        auto [elmat_Tc, elmat_Tt] = elmat_T.SplitCols (ndof_conforming);

        // T_c solves A @ T_c = B,
        elmat_Tc = elmat_A_inv * elmat_B;

        if (fes_ip)
          {
            Matrix<SCAL> Tc (fes_ip_sqinv * elmat_Tc);
            elmat_Tc = Tc;
            FlatMatrix<SCAL, ColMajor> Vtemp (V.Width (), V.Height (), lh);
            Vtemp = fes_ip_sqinv * V;
            V.Assign (Vtemp);
          }
        elmat_Tt = V.Cols (ndof - ndof_trefftz_i, ndof);

        // Save the particular-solution map in the same local ordering.
        {
          auto elmat_T_inv = elmat_A_inv.Cols (ndof_conforming,
                                               ndof_conforming + ndof_test);

          if (fes_ip)
            {
              FlatMatrix<SCAL> elmat_T_inv_phys (ndof, ndof_test, lh);
              elmat_T_inv_phys = fes_ip_sqinv * elmat_T_inv;
              elmat_T_inv = elmat_T_inv_phys;
            }

          work.inverses[nr] = elmat_T_inv;
        }

        if (stats)
          {
            const lock_guard<mutex> lock (work.stats_mutex);
            auto diag = elmat_A.Diag (0);
            if (work.sing_val_avg.Size () < diag.Size ())
              {
                const size_t old_size = work.sing_val_avg.Size ();
                Vector<double> avg (diag.Size ()), maxval (diag.Size ()),
                    minval (diag.Size ()), count (diag.Size ());
                avg = 0.;
                maxval = 0.;
                minval = DBL_MAX;
                count = 0.;
                avg.Range (0, old_size) = work.sing_val_avg;
                maxval.Range (0, old_size) = work.sing_val_max;
                minval.Range (0, old_size) = work.sing_val_min;
                count.Range (0, old_size) = work.sing_val_count;
                work.sing_val_avg = std::move (avg);
                work.sing_val_max = std::move (maxval);
                work.sing_val_min = std::move (minval);
                work.sing_val_count = std::move (count);
              }
            for (size_t i = 0; i < diag.Size (); i++)
              {
                work.sing_val_count[i] += 1.;
                work.sing_val_avg[i] += abs (diag (i));
                work.sing_val_max[i]
                    = max (work.sing_val_max[i], abs (diag (i)));
                work.sing_val_min[i]
                    = min (work.sing_val_min[i], abs (diag (i)));
              }
            Matrix<SCAL> test = elmat_A_copy * elmat_Tc - elmat_B;
            for (auto &cr : test.AsVector ())
              if (abs (cr) > work.conformity_residual)
                work.conformity_residual = abs (cr);
            test = elmat_A_copy * elmat_Tt;
            for (auto &tr : test.AsVector ())
              if (abs (tr) > work.trefftz_residual)
                work.trefftz_residual = abs (tr);
          }
      }

    if (trhs)
      {
        FlatVector<SCAL> rhs (ndof_test, lh), sol (ndof, lh);
        rhs = SCAL (0.0);
        add_vector (rhs, dofs_test, lh);
        sol = *work.inverses[nr] * rhs;
        if (work.patchwise)
          work.particular_solution->AddIndirect (dofs, sol, true);
        else
          work.particular_solution->SetIndirect (dofs, sol);
      }

    if (work.patchwise)
      {
        data.active = true;
        data.ndof_trefftz = ndof_trefftz_i;
        data.shared = local_data[nr].shared;
        local_data[nr] = std::move (data);
      }
    else
      {
        local_ndofs_trefftz[nr] = ndof_trefftz_i;
        // Restore hidden rows and ignored identity columns only when needed.
        if (work.any_fes_has_inactive_dofs || ignoredofs)
          {
            Array<DofId> all_dofs;
            get_dofs (*fes, all_dofs);
            size_t nignored = 0;
            for (DofId d : all_dofs)
              if (IsIgnoredRegularDof (ignoredofs, d))
                nignored++;
            if (ndof != all_dofs.Size () || nignored)
              {
                Matrix<SCAL> full (all_dofs.Size (),
                                   elmat_T.Width () + nignored);
                full = SCAL (0.0);
                size_t row = 0, col = elmat_T.Width ();
                for (size_t i : Range (all_dofs))
                  if (IsIgnoredRegularDof (ignoredofs, all_dofs[i]))
                    full (i, col++) = SCAL (1.0);
                  else if (!IsInactiveOrHiddenDof (*fes, all_dofs[i]))
                    full.Cols (elmat_T.Width ()).Row (i) = elmat_T.Row (row++);
                elmat_T = std::move (full);
              }
          }
      }
    work.local_matrices[nr] = std::move (elmat_T);
  }

  template <typename SCAL> void TrefftzEmbedding::EmbTrefftz ()
  {

    static Timer t ("EmbTrefftz: EmbTrefftz");
    RegionTimer reg (t);
    EmbeddingWork<SCAL> work;
    static_assert (std::is_same_v<double, SCAL>
                       || std::is_same_v<Complex, SCAL>,
                   "SCAL must be double or Complex");

    vb = ma->GetNE (VOL) > 0 ? VOL : BND;
    const size_t num_elements = ma->GetNE (vb);
    // The space supplies its partition; ordinary spaces stay elementwise.
#ifdef NGSTREFFTZ_USE_XFEM
    Table<size_t> patch_facets;
    if (patches)
      patch_facets = GetInnerPatchFacets (*ma, patches->elements);
#endif
    work.patchwise = bool (patches);
    size_t num_local
        = work.patchwise ? patches->elements.Size () : num_elements;
    if (work.patchwise)
      {
        local_data.SetSize (num_local);
        for (size_t p : Range (num_local))
          local_data[p].shared = patches->elements[p].Size () > 1;
      }
    local_ndofs_trefftz.SetSize (num_elements);
    local_ndofs_trefftz = 0;

    // #TODO what is a good size for the local heap?
    // For the moment: large enough constant size.
    LocalHeap clh = LocalHeap (100 * 1000 * 1000, "embt", true);

    // calculate the integrators for the three bilinear forms,
    // each for VOL, BND, BBND, BBBND, hence 4 arrays per bilnear form
    if (top)
      calculateBilinearFormIntegrators (*top, work.op_integrators, ma);
    if (cop && crhs)
      {
        calculateBilinearFormIntegrators (*cop, work.cop_lhs_integrators, ma);
        calculateBilinearFormIntegrators (*crhs, work.cop_rhs_integrators, ma);
      }
    if (fes_ip)
      calculateBilinearFormIntegrators (*fes_ip, work.fes_ip_integrators, ma);

    work.local_matrices.SetSize (num_local);
    work.inverses.SetSize (num_local);
    work.particular_solution = make_shared<VVector<SCAL>> (fes->GetNDof ());
    *work.particular_solution = SCAL (0.0);
    if (trhs)
      CalculateLinearFormIntegrators (*trhs, work.lfis, ma, bool (patches));

    work.test_has_inactive_dofs = fesHasInactiveDofs (*fes_test);
    work.any_fes_has_inactive_dofs
        = fesHasInactiveDofs (*fes) || work.test_has_inactive_dofs
          || ((fes_conformity) ? fesHasInactiveDofs (*fes_conformity) : false);

    auto defined_on = [&] (const Ngs_Element &el) {
      return (top && bfIsDefinedOnElement (*top, el, ma))
             || (cop && bfIsDefinedOnElement (*cop, el, ma))
             || (crhs && bfIsDefinedOnElement (*crhs, el, ma));
    };

    if (!work.patchwise)
      ma->IterateElements (vb, clh, [&] (Ngs_Element el, LocalHeap &lh) {
        ElementId ei (el);
        auto get_dofs = [&] (const FESpace &space, Array<DofId> &dofs) {
          space.GetDofNrs (ei, dofs);
        };
        auto add_matrix
            = [&] (auto mat, const auto &integrators, const FESpace &trial,
                   const FESpace &test, FlatArray<DofId>, FlatArray<DofId>,
                   LocalHeap &lh) {
                AddIntegrationToElementMatrix (mat, integrators, *ma, ei,
                                               trial, test, lh);
              };
        auto add_vector = [&] (FlatVector<SCAL> vec, FlatArray<DofId>,
                               LocalHeap &lh) {
          vec += AssembleElementRhs<SCAL> (work.lfis, *fes_test, ei, *ma,
                                           work.test_has_inactive_dofs, lh);
        };
        if (defined_on (el))
          ComputeLocalEmbedding<SCAL> (el.Nr (), work, get_dofs, add_matrix,
                                       add_vector, lh);
      });
    else
      {
#ifdef NGSTREFFTZ_USE_XFEM
        ParallelForRange (Range (num_local), [&] (IntRange range) {
          LocalHeap lh = clh.Split ();
          for (size_t nr : range)
            {
              HeapReset reset (lh);
              PatchAssembly local{ *ma, patches->elements[nr],
                                   patch_facets[nr] };
              bool defined = false;
              for (size_t el : local.elements)
                if (defined_on (ma->GetElement (ElementId (vb, el))))
                  {
                    defined = true;
                    break;
                  }
              if (!defined)
                continue;
              auto get_dofs = [&] (const FESpace &space, Array<DofId> &dofs) {
                local.GetDofNrs (space, dofs);
              };
              auto add_matrix =
                  [&] (auto mat, const auto &integrators, const FESpace &trial,
                       const FESpace &test, FlatArray<DofId> dofs,
                       FlatArray<DofId> test_dofs, LocalHeap &lh) {
                    local.AddMatrix<SCAL> (mat, integrators, trial, test, dofs,
                                           test_dofs, lh);
                  };
              auto add_vector
                  = [&] (FlatVector<SCAL> vec, FlatArray<DofId> test_dofs,
                         LocalHeap &lh) {
                      local.AddVector (vec, work.lfis[VOL], *fes_test,
                                       test_dofs, lh);
                    };
              ComputeLocalEmbedding<SCAL> (nr, work, get_dofs, add_matrix,
                                           add_vector, lh);
            }
        });
#else
        throw Exception ("Patch embedding requires xfem support");
#endif
      }
    size_t offset = 0;
    for (auto &data : local_data)
      {
        data.trefftz_dof_offset = offset;
        offset += data.ndof_trefftz;
      }

    if (stats)
      {
        for (size_t i : Range (work.sing_val_avg.Size ()))
          work.sing_val_avg[i] /= work.sing_val_count[i];
        (*stats)["singavg"] = Vector<double> (work.sing_val_avg);
        (*stats)["singmax"] = Vector<double> (work.sing_val_max);
        (*stats)["singmin"] = Vector<double> (work.sing_val_min);
        (*stats)["conformity_residual"]
            = Vector<double> ({ work.conformity_residual });
        (*stats)["trefftz_residual"]
            = Vector<double> ({ work.trefftz_residual });
      }

    if (work.patchwise)
      {
        Array<optional<Matrix<SCAL>>> element_matrices;
        AssemblePatchEmbeddings (
            work.local_matrices, local_data, *patches, *fes, fes_conformity,
            ignoredofs, element_matrices, tdof_nrs, local_ndofs_trefftz);
        work.local_matrices = std::move (element_matrices);
      }
    else
      tdof_nrs
          = ElementEmbeddingDofs (work.local_matrices, local_ndofs_trefftz,
                                  *fes, fes_conformity, ignoredofs);

    if constexpr (std::is_same_v<double, SCAL>)
      {
        etmats = std::move (work.local_matrices);
        local_inverses = std::move (work.inverses);
      }
    else
      {
        etmatsc = std::move (work.local_matrices);
        local_inverses_c = std::move (work.inverses);
      }
    psol = std::move (work.particular_solution);
  }

  TrefftzEmbedding::TrefftzEmbedding (
      shared_ptr<SumOfIntegrals> _top, shared_ptr<SumOfIntegrals> _trhs,
      shared_ptr<SumOfIntegrals> _cop, shared_ptr<SumOfIntegrals> _crhs,
      size_t _ndof_trefftz, double _eps, shared_ptr<FESpace> _fes,
      shared_ptr<FESpace> _fes_test, shared_ptr<FESpace> _fes_conformity,
      shared_ptr<SumOfIntegrals> _fes_ip, shared_ptr<BitArray> _ignoredofs,
      shared_ptr<std::map<std::string, Vector<double>>> _stats)
      : top (_top), trhs (_trhs), cop (_cop), crhs (_crhs), fes_ip (_fes_ip),
        ignoredofs (_ignoredofs), stats (_stats)
  {
    if (_ndof_trefftz == std::numeric_limits<size_t>::max ())
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
    patches = GetDofPatches (*fes);
    if (fes->IsComplex ())
      EmbTrefftz<Complex> ();
    else
      EmbTrefftz<double> ();
  }

  size_t TrefftzEmbedding::GetNDofTrefftz () const noexcept
  {
    size_t ndof = 0;
    if (patches)
      for (const auto &data : local_data)
        ndof += data.ndof_trefftz;
    else
      for (size_t n : local_ndofs_trefftz)
        ndof += n;
    return ndof;
  }

  template <typename SCAL>
  shared_ptr<const BaseVector> TrefftzEmbedding::ComputeParticularSolution (
      shared_ptr<const BaseVector> rhsvec,
      const Array<optional<Matrix<SCAL>>> &inverses,
      const SumOfIntegrals *rhsform) const
  {
    auto solution = make_shared<VVector<SCAL>> (fes->GetNDof ());
    *solution = SCAL (0.0);
    if (!rhsvec && !rhsform)
      return solution;
    LocalHeap lh (100 * 1000 * 1000, "embt_psol", true);
    Array<shared_ptr<LinearFormIntegrator>> lfis[4];
    if (rhsform)
      CalculateLinearFormIntegrators (*rhsform, lfis, ma, bool (patches));
    const bool test_has_inactive_dofs
        = rhsform && fesHasInactiveDofs (*fes_test);
    auto apply = [&] (size_t nr, FlatArray<DofId> dofs,
                      FlatArray<DofId> test_dofs, const auto &add_vector,
                      LocalHeap &lh) {
      FlatVector<SCAL> rhs (test_dofs.Size (), lh), sol (dofs.Size (), lh);
      rhs = SCAL (0.0);
      if (rhsvec)
        rhsvec->GetIndirect (test_dofs, rhs);
      else
        add_vector (rhs, test_dofs, lh);
      sol = *inverses[nr] * rhs;
      if (patches)
        solution->AddIndirect (dofs, sol, true);
      else
        solution->SetIndirect (dofs, sol);
    };
    if (!patches)
      ma->IterateElements (vb, lh, [&] (Ngs_Element el, LocalHeap &lh) {
        ElementId ei (el);
        if (!inverses[ei.Nr ()])
          return;
        Array<DofId> dofs, test_dofs;
        fes->GetDofNrs (ei, dofs, VISIBLE_DOF);
        fes_test->GetDofNrs (ei, test_dofs, VISIBLE_DOF);
        size_t n = 0;
        for (DofId d : dofs)
          if (!IsInactiveOrHiddenDof (*fes, d)
              && !IsIgnoredRegularDof (ignoredofs, d))
            dofs[n++] = d;
        dofs.SetSize (n);
        n = 0;
        for (DofId d : test_dofs)
          if (!IsInactiveOrHiddenDof (*fes_test, d))
            test_dofs[n++] = d;
        test_dofs.SetSize (n);
        auto add_vector
            = [&] (FlatVector<SCAL> rhs, FlatArray<DofId>, LocalHeap &lh) {
                rhs += AssembleElementRhs<SCAL> (lfis, *fes_test, ei, *ma,
                                                 test_has_inactive_dofs, lh);
              };
        apply (ei.Nr (), dofs, test_dofs, add_vector, lh);
      });
    else
      {
#ifdef NGSTREFFTZ_USE_XFEM
        ParallelForRange (Range (local_data.Size ()), [&] (IntRange range) {
          LocalHeap local_heap = lh.Split ();
          for (size_t nr : range)
            if (local_data[nr].active)
              {
                HeapReset reset (local_heap);
                const auto &data = local_data[nr];
                PatchAssembly local{ *ma, patches->elements[nr], {} };
                auto add_vector = [&] (FlatVector<SCAL> rhs,
                                       FlatArray<DofId> test_dofs,
                                       LocalHeap &lh) {
                  local.AddVector (rhs, lfis[VOL], *fes_test, test_dofs, lh);
                };
                apply (nr, data.trial_dofs, data.test_dofs, add_vector,
                       local_heap);
              }
        });
#endif
      }
    return solution;
  }

  shared_ptr<const BaseVector> TrefftzEmbedding::GetParticularSolution () const
  {
    return psol;
  }

  shared_ptr<const BaseVector> TrefftzEmbedding::GetParticularSolution (
      shared_ptr<SumOfIntegrals> _trhs) const
  {
    if (fes->IsComplex ())
      return ComputeParticularSolution<Complex> (nullptr, local_inverses_c,
                                                 _trhs.get ());
    return ComputeParticularSolution<double> (nullptr, local_inverses,
                                              _trhs.get ());
  }

  shared_ptr<const BaseVector> TrefftzEmbedding::GetParticularSolution (
      shared_ptr<const BaseVector> _trhsvec) const
  {
    if (fes->IsComplex ())
      return ComputeParticularSolution (_trhsvec, local_inverses_c);
    return ComputeParticularSolution (_trhsvec, local_inverses);
  }

  shared_ptr<const BaseMatrix> TrefftzEmbedding::GetEmbedding () const
  {
    size_t ndof = GetNDofTrefftz ();
    if (fes_conformity)
      ndof += fes_conformity->GetNDof ();
    if (ignoredofs)
      ndof += ignoredofs->NumSet ();
    if (fes->IsComplex ())
      return Elmats2Sparse (etmatsc, tdof_nrs, *fes, ndof, ignoredofs,
                            bool (patches));
    return Elmats2Sparse (etmats, tdof_nrs, *fes, ndof, ignoredofs,
                          bool (patches));
  }

  template <typename SCAL>
  shared_ptr<BaseVector>
  EmbedElementMatrices (const BaseVector &tvec,
                        const Array<optional<Matrix<SCAL>>> &etmats,
                        const Table<DofId> &tdof_nrs, const FESpace &fes,
                        VorB vb, bool add_contributions,
                        shared_ptr<const BitArray> ignoredofs)
  {
    auto vec = make_shared<VVector<SCAL>> (fes.GetNDof ());
    *vec = SCAL (0.0);
    LocalHeap lh (sizeof (Complex) * 10000, "embt", true);
    // Base-space coloring also covers globally shared ignored DOFs.
    IterateElements (fes, vb, lh, [&] (auto ei, LocalHeap &lh) {
      if (!etmats[ei.Nr ()])
        return;
      Array<DofId> dofs;
      fes.GetDofNrs (ei, dofs);
      auto tdofs = tdof_nrs[ei.Nr ()];
      FlatVector<SCAL> telvec (tdofs.Size (), lh), elvec (dofs.Size (), lh);
      tvec.GetIndirect (tdofs, telvec);
      elvec = *etmats[ei.Nr ()] * telvec;
      if (add_contributions)
        vec->AddIndirect (dofs, elvec);
      else
        vec->SetIndirect (dofs, elvec);
    });
    // An ignored DOF is one global identity column, even if shared.
    if (ignoredofs)
      {
        size_t next = tvec.Size () - ignoredofs->NumSet ();
        for (size_t d : Range (ignoredofs->Size ()))
          if (ignoredofs->Test (d))
            vec->template FV<SCAL> ()[d] = tvec.FV<SCAL> ()[next++];
      }
    return vec;
  }

  shared_ptr<BaseVector>
  TrefftzEmbedding::Embed (const shared_ptr<const BaseVector> tvec) const
  {
    bool add_contributions = bool (fes_conformity) && !patches;
    if (fes->IsComplex ())
      return EmbedElementMatrices (*tvec, etmatsc, tdof_nrs, *fes, vb,
                                   add_contributions, ignoredofs);
    return EmbedElementMatrices (*tvec, etmats, tdof_nrs, *fes, vb,
                                 add_contributions, ignoredofs);
  }

  shared_ptr<GridFunction>
  TrefftzEmbedding::Embed (const shared_ptr<const GridFunction> tgfu) const
  {
    const shared_ptr<GridFunction> gfu
        = CreateGridFunction (this->fes, tgfu->GetName (), Flags ());
    gfu->Update ();
    *gfu->GetVectorPtr () = *this->Embed (tgfu->GetVectorPtr ());
    return gfu;
  }

  ////////////////////////// EmbeddedTrefftzFES ///////////////////////////

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

  EmbeddedTrefftzFES::EmbeddedTrefftzFES (shared_ptr<TrefftzEmbedding> aemb)
      : FESpace (aemb->GetFES ()->GetMeshAccess (),
                 aemb->GetFES ()->GetFlags (), false),
        emb (aemb), basefes (aemb->GetFES ()), etmats (emb->GetEtmats ()),
        etmatsc (emb->GetEtmatsC ()), fes_conformity (emb->GetConfSpace ()),
        ignoredofs (emb->GetIgnoredDofs ())
  {
    this->name = "EmbeddedTrefftzFES(" + basefes->GetClassName () + ")";
    this->type = "embt";
    this->needs_transform_vec = true;
    this->iscomplex = basefes->IsComplex ();

    for (auto vb : { VOL, BND, BBND, BBBND })
      {
        this->evaluator[vb] = basefes->GetEvaluator (vb);
        this->flux_evaluator[vb] = basefes->GetFluxEvaluator (vb);
      }

    if (this->IsComplex ())
      etmatsc_inv.SetSize (this->ma->GetNE (emb->VB ()));
    else
      etmats_inv.SetSize (this->ma->GetNE (emb->VB ()));

    this->Update ();

    if (this->order_policy == ORDER_POLICY::VARIABLE_ORDER)
      {
        this->GetMeshAccess ()->IterateElements (VOL, [&] (ElementId ei) {
          NodeId ni (NODE_TYPE::NT_ELEMENT, ei.Nr ());
          this->SetOrder (ni, basefes->GetOrder (ni));
        });
        this->Update ();
      }

    this->FinalizeUpdate ();

    if (fes_conformity)
      copyBitArray (this->GetFreeDofs (), fes_conformity->GetFreeDofs ());
  }

  void EmbeddedTrefftzFES::Update ()
  {
    basefes->Update ();
    if (fes_conformity)
      fes_conformity->Update ();

    if (this->order_policy == ORDER_POLICY::VARIABLE_ORDER)
      {
        this->GetMeshAccess ()->IterateElements (VOL, [&] (ElementId ei) {
          NodeId ni (NODE_TYPE::NT_ELEMENT, ei.Nr ());
          this->SetOrder (ni, basefes->GetOrder (ni));
        });
      }

    for (auto vb : { VOL, BND, BBND, BBBND })
      {
        this->evaluator[vb] = basefes->GetEvaluator (vb);
        this->flux_evaluator[vb] = basefes->GetFluxEvaluator (vb);
      }
    this->iscomplex = basefes->IsComplex ();

    this->UpdateDofTables ();
    this->UpdateCouplingDofArray ();
    this->UpdateFreeDofs ();
    FESpace::Update ();
  }

  void EmbeddedTrefftzFES::UpdateDofTables ()
  {
    size_t ndof_conforming = fes_conformity ? fes_conformity->GetNDof () : 0;
    size_t ignored_dofs = ignoredofs ? ignoredofs->NumSet () : 0;
    const size_t ndof_trefftz = emb->GetNDofTrefftz ();
    this->SetNDof (ndof_conforming + ignored_dofs + ndof_trefftz);

    if (this->IsComplex ())
      etmatsc_inv.SetSize (this->ma->GetNE (emb->VB ()));
    else
      etmats_inv.SetSize (this->ma->GetNE (emb->VB ()));
  }

  void EmbeddedTrefftzFES::UpdateFreeDofs ()
  {
    this->free_dofs = make_shared<BitArray> (this->GetNDof ());
    this->free_dofs->Set ();

    size_t ndof_conforming = fes_conformity ? fes_conformity->GetNDof () : 0;
    if (fes_conformity)
      {
        auto conf_free = fes_conformity->GetFreeDofs ();
        for (size_t i = 0;
             i < min<size_t> (ndof_conforming, conf_free->Size ()); i++)
          if (!conf_free->Test (i))
            this->free_dofs->Clear (i);
      }

    if (ignoredofs)
      {
        const size_t ndof_trefftz = emb->GetNDofTrefftz ();

        auto base_free = basefes->GetFreeDofs ();
        size_t idof = ndof_conforming + ndof_trefftz;
        for (size_t i = 0; i < ignoredofs->Size (); i++)
          if (ignoredofs->Test (i))
            {
              if (i < base_free->Size () && !base_free->Test (i))
                this->free_dofs->Clear (idof);
              idof++;
            }
      }
  }

  void EmbeddedTrefftzFES::UpdateCouplingDofArray ()
  {
    const size_t ndof_conformity
        = (fes_conformity) ? fes_conformity->GetNDof () : 0;

    const size_t ndof_trefftz = emb->GetNDofTrefftz ();

    size_t ignored_dofs = 0;
    if (ignoredofs)
      ignored_dofs = ignoredofs->NumSet ();

    // The conformity dofs might overlap.
    // Overall, they add up to dofs in the conformity space.
    const size_t new_ndof = ignored_dofs + ndof_conformity + ndof_trefftz;
    this->SetNDof (new_ndof);
    this->ctofdof.SetSize (new_ndof);
    // Global dof ordering (consistent with AssemblePatchEmbeddings):
    //   [0 .. ndof_conformity-1]                                  conforming
    //   dofs [ndof_conformity .. ndof_conformity+ndof_trefftz-1] element-wise
    //   Trefftz dofs [ndof_conformity+ndof_trefftz .. end] ignored base dofs
    // 1) Conformity block (preserve base coupling types)
    if (fes_conformity)
      for (size_t i = 0; i < ndof_conformity; i++)
        this->ctofdof[i] = fes_conformity->GetDofCouplingType (i);
    // 2) Trefftz block (ordinary element modes are local).
    this->ctofdof.Range (ndof_conformity, ndof_conformity + ndof_trefftz)
        = LOCAL_DOF;
    for (const auto &data : emb->GetLocalData ())
      for (size_t j : Range (data.ndof_trefftz))
        this->ctofdof[ndof_conformity + data.trefftz_dof_offset + j]
            = data.shared ? WIREBASKET_DOF : LOCAL_DOF;
    // 3) Ignored block (preserve base coupling types)
    if (ignoredofs)
      {
        size_t idof = ndof_conformity + ndof_trefftz;
        for (size_t i = 0; i < ignoredofs->Size (); i++)
          if (ignoredofs->Test (i))
            this->ctofdof[idof++] = basefes->GetDofCouplingType (i);
      }
  }

  void EmbeddedTrefftzFES::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize0 ();
    if (!basefes->DefinedOn (ei) || ei.VB () != emb->VB ())
      return;
    const FlatArray<DofId> tdofnrs = this->emb->GetTDofNrs (ei.Nr ());

    basefes->GetDofNrs (ei, dnums);
    const size_t base_ndof = dnums.Size ();
    if (tdofnrs.Size () > base_ndof)
      {
        dnums.SetSize (tdofnrs.Size ());
        for (size_t i = base_ndof; i < tdofnrs.Size (); i++)
          dnums[i] = tdofnrs[i];
      }
    for (size_t i = 0; i < dnums.Size (); i++)
      dnums[i] = i < tdofnrs.Size () ? tdofnrs[i] : NO_DOF_NR_CONDENSE;
  }

  FiniteElement &
  EmbeddedTrefftzFES::GetFE (ElementId ei, Allocator &alloc) const
  {
    return basefes->GetFE (ei, alloc);
  }

  optional<FlatMatrix<double>>
  EmbeddedTrefftzFES::GetEtmatInv (size_t idx) const
  {
    std::call_once (this->etmats_inv_computed, [&] () {
      this->GetMeshAccess ()->IterateElements (emb->VB (), [&] (ElementId ei) {
        const optional<Matrix<double>> &etmat = emb->GetEtmat (ei.Nr ());
        this->etmats_inv[ei.Nr ()]
            = (etmat) ? make_optional (getPseudoInverse (*etmat, 0)) : nullopt;
      });
    });
    return this->etmats_inv[idx];
  }

  optional<FlatMatrix<Complex>>
  EmbeddedTrefftzFES::GetEtmatCInv (size_t idx) const
  {
    std::call_once (this->etmats_inv_computed, [&] () {
      this->GetMeshAccess ()->IterateElements (emb->VB (), [&] (ElementId ei) {
        const optional<Matrix<Complex>> &etmat = emb->GetEtmatC (ei.Nr ());
        this->etmatsc_inv[ei.Nr ()]
            = (etmat) ? make_optional (getPseudoInverse (*etmat, 0)) : nullopt;
      });
    });
    return this->etmatsc_inv[idx];
  }

  ProxyNode EmbeddedTrefftzFES::MakeProxyFunction (
      bool testfunction,
      const function<shared_ptr<ProxyFunction> (shared_ptr<ProxyFunction>)>
          &addblock) const
  {
    ProxyNode pn = basefes->MakeProxyFunction (testfunction, addblock);
    auto self = dynamic_pointer_cast<FESpace> (
        const_cast<EmbeddedTrefftzFES *> (this)->shared_from_this ());
    pn.SetFESpace (self);
    return pn;
  }

  template <typename SCAL>
  void T_VTransformM (SliceMatrix<SCAL> mat, const TRANSFORM_TYPE type,
                      const Matrix<SCAL> &elmat)
  {
    static Timer timer ("EmbTrefftz: TransformM");
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

  void EmbeddedTrefftzFES::VTransformMR (ElementId ei, SliceMatrix<double> mat,
                                         TRANSFORM_TYPE type) const
  {
    if (!etmats[ei.Nr ()])
      {
        mat = 0.0;
        return;
      }
    T_VTransformM (mat, type, *etmats[ei.Nr ()]);
  }

  void
  EmbeddedTrefftzFES::VTransformMC (ElementId ei, SliceMatrix<Complex> mat,
                                    TRANSFORM_TYPE type) const
  {
    if (!etmatsc[ei.Nr ()])
      {
        mat = 0.0;
        return;
      }
    T_VTransformM (mat, type, *etmatsc[ei.Nr ()]);
  }

  template <typename SCAL>
  void T_VTransformV (SliceVector<SCAL> &vec, const TRANSFORM_TYPE type,
                      const Matrix<SCAL> &elmat,
                      optional<FlatMatrix<SCAL>> elmat_inv)
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
        Vector<SCAL> tmp_vec ((*elmat_inv).Height ());
        tmp_vec = (*elmat_inv) * vec;
        vec = tmp_vec;
      }
    else
      {
        stringstream err;
        err << "VTransformV is not implemented for TRANSFORM_TYPE " << type;
        throw std::invalid_argument (err.str ());
      }
  }

  void EmbeddedTrefftzFES::VTransformVR (ElementId ei, SliceVector<double> vec,
                                         TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: TransformV");
    RegionTimer reg (timer);

    if (!emb->GetEtmat (ei.Nr ()))
      {
        vec = 0.0;
        return;
      }

    T_VTransformV (vec, type, *(emb->GetEtmat (ei.Nr ())),
                   (type == TRANSFORM_SOL_INVERSE)
                       ? this->GetEtmatInv (ei.Nr ())
                       : nullopt);
  }

  void
  EmbeddedTrefftzFES::VTransformVC (ElementId ei, SliceVector<Complex> vec,
                                    TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: TransformV");
    RegionTimer reg (timer);

    if (!emb->GetEtmatC (ei.Nr ()))
      {
        vec = 0.0;
        return;
      }

    T_VTransformV (vec, type, *(emb->GetEtmatC (ei.Nr ())),
                   (type == TRANSFORM_SOL_INVERSE)
                       ? this->GetEtmatCInv (ei.Nr ())
                       : nullopt);
  }

}

string EmbeddedTrefftzFES::GetClassName () const { return this->name; }

////////////////////////// python interface ///////////////////////////

#ifdef NGS_PYTHON
void ExportEmbTrefftz (py::module m)
{
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

                Spaces may supply independent trial-DOF groups and their element supports.
                The same DOFs select test rows when fes_test is fes; otherwise test rows
                are collected over the patch elements. ignoredofs only removes trial columns.
            )mydelimiter")
      .def (
          // could be py::init<..,optional<py::dict>>() but we translate to map
          py::init ([] (shared_ptr<SumOfIntegrals> top,
                        shared_ptr<SumOfIntegrals> trhs,
                        shared_ptr<SumOfIntegrals> cop,
                        shared_ptr<SumOfIntegrals> crhs, size_t ndof_trefftz,
                        double eps, shared_ptr<FESpace> fes,
                        shared_ptr<FESpace> fes_test,
                        shared_ptr<FESpace> fes_conformity,
                        shared_ptr<SumOfIntegrals> fes_ip,
                        shared_ptr<BitArray> ignoredofs,
                        optional<py::dict> pystats) {
            shared_ptr<std::map<std::string, Vector<double>>> stats = nullptr;
            if (pystats)
              stats = make_shared<std::map<std::string, Vector<double>>> ();
            std::shared_ptr<TrefftzEmbedding> emb
                = std::make_shared<TrefftzEmbedding> (
                    top, trhs, cop, crhs, ndof_trefftz, eps, fes, fes_test,
                    fes_conformity, fes_ip, ignoredofs, stats);
            if (pystats)
              for (auto const &x : *stats)
                (*pystats)[py::cast (x.first)] = py::cast (x.second);
            return emb;
          }),
          R"mydelimiter(
                Constructs a new Trefftz embedding object.

                 :param top: the differential operation. Can be None
                 :param trhs: right hand side of the var. formulation
                 :param cop: left hand side of the conformity operation
                 :param crhs: right hand side of the conformity operation
                 :param eps: cutoff for singular values from the SVD of the local operator.
                        values below eps are considered zero and therefore in the kernel of `top`.
                        (default: 0.0)
                 :param ndof_trefftz: fixes the number of degrees of freedom per element or patch
                     that are to be considered in the Trefftz space generated by `top`
                     (i.e. the local dimension of the kernel of `top` on one element or patch)
                     cannot be used together with `eps` (default: 0)
                 :param fes: the finite element space of `top` (optional, determined
                     from `top` if not given)
                 :param fes_test: the finite element test space of `top` (optional,
                     determined from `top` if not given)
                 :param fes_conformity: finite element space of the conformity operation (optional,
                     determined from `cop` if not given)
                 :param ignoredofs: BitArray of dofs from fes to be ignored in the embedding
                 :param stats: optional dictionary to store statistics about the singular values,
                     input dictionary is modified
            )mydelimiter",
          py::arg ("top") = nullptr, py::arg ("trhs") = nullptr,
          py::arg ("cop") = nullptr, py::arg ("crhs") = nullptr,
          py::arg ("ndof_trefftz") = std::numeric_limits<size_t>::max (),
          py::arg ("eps") = 0.0, py::arg ("fes") = nullptr,
          py::arg ("fes_test") = nullptr, py::arg ("fes_conformity") = nullptr,
          py::arg ("fes_ip") = nullptr, py::arg ("ignoredofs") = nullptr,
          py::arg ("stats") = nullopt)
      .def ("Embed",
            static_cast<shared_ptr<BaseVector> (ngcomp::TrefftzEmbedding::*) (
                const shared_ptr<const BaseVector>) const> (
                &ngcomp::TrefftzEmbedding::Embed),
            "Embed a Trefftz GridFunction Vector into the underlying FESpace")
      .def (
          "Embed",
          static_cast<shared_ptr<GridFunction> (ngcomp::TrefftzEmbedding::*) (
              const shared_ptr<const GridFunction>) const> (
              &ngcomp::TrefftzEmbedding::Embed),
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
            "FESpace, given a trhs")
      .def ("GetParticularSolution",
            static_cast<shared_ptr<const BaseVector> (TrefftzEmbedding::*) (
                shared_ptr<const BaseVector>) const> (
                &ngcomp::TrefftzEmbedding::GetParticularSolution),
            "Particular solution as GridFunction vector of the underlying "
            "FESpace, given a trhs as vector")
      .def_property_readonly (
          "fes",
          [] (shared_ptr<TrefftzEmbedding> emb) { return emb->GetFES (); })
      .def_property_readonly ("fes_test",
                              [] (shared_ptr<TrefftzEmbedding> emb) {
                                return emb->GetTestSpace ();
                              })
      .def_property_readonly ("fes_conformity",
                              [] (shared_ptr<TrefftzEmbedding> emb) {
                                return emb->GetConfSpace ();
                              });

  auto pyspace = ngcomp::ExportFESpace<ngcomp::EmbeddedTrefftzFES> (
      m, "EmbeddedTrefftzFES");

  pyspace.def (py::init ([] (shared_ptr<ngcomp::TrefftzEmbedding> emb) {
                 auto nfes = make_shared<ngcomp::EmbeddedTrefftzFES> (emb);
                 nfes->Update ();
                 connect_auto_update (nfes.get ());
                 return nfes;
               }),
               R"mydelimiter(
                Constructs a new EmbeddedTrefftzFESpace from a TrefftzEmbedding.

                :param emb: TrefftzEmbedding object that defines the embedding
               )mydelimiter",
               py::arg ("emb"));

  pyspace.def ("GetEmbedding", &ngcomp::EmbeddedTrefftzFES::GetEmbedding,
               "Get the TrefftzEmbedding");
  pyspace.def_property_readonly ("emb", &EmbeddedTrefftzFES::GetEmbedding);
  pyspace.def_property_readonly ("base", &EmbeddedTrefftzFES::GetBaseFESpace);
  pyspace.def_property_readonly (
      "components", [] (shared_ptr<EmbeddedTrefftzFES> self) {
        py::list lst;
        auto base = self->GetBaseFESpace ();
        if (auto comp = dynamic_pointer_cast<CompoundFESpace> (base))
          for (auto sp : comp->Spaces ())
            lst.append (sp);
        return lst;
      });
}

#endif // NGS_PYTHON
