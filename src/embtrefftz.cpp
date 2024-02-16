#include "embtrefftz.hpp"
#include "monomialfespace.hpp"
#include <cfloat>

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

  template <class SCAL>
  void GetSVD (SliceMatrix<SCAL> A, SliceMatrix<SCAL, ColMajor> U,
               SliceMatrix<SCAL, ColMajor> V)
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
    A = 0.0;
    // A.Diag(0)=AA.Diag();
    for (int i = 0; i < min (A.Width (), A.Height ()); i++)
      A (i, i) = AA (i, i);
  }

  template void
  GetSVD<double> (SliceMatrix<double> A, SliceMatrix<double, ColMajor> U,
                  SliceMatrix<double, ColMajor> V);

  template <>
  void
  GetSVD<Complex> (SliceMatrix<Complex> A, SliceMatrix<Complex, ColMajor> U,
                   SliceMatrix<Complex, ColMajor> V)
  {
    Matrix<Complex, ColMajor> AA = A;
#ifdef NGSTREFFTZ_USE_LAPACK
    LapackSVD (AA, U, V);
#else
    throw Exception ("Need Lapack for complex SVD");
#endif
    A = 0.0;
    for (int i = 0; i < min (A.Width (), A.Height ()); i++)
      A (i, i) = AA (i, i);
  }
}

namespace ngcomp
{

  template <class SCAL>
  std::tuple<vector<shared_ptr<Matrix<SCAL>>>, shared_ptr<BaseVector>>
  EmbTrefftz (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
              shared_ptr<SumOfIntegrals> lf, double eps,
              shared_ptr<FESpace> test_fes, int tndof, bool getrange,
              std::map<std::string, Vector<SCAL>> *stats)
  {
    static Timer svdtt ("svdtrefftz");
    RegionTimer reg (svdtt);
    LocalHeap lh (1000 * 1000 * 1000);

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
    for (auto icf : bf->icfs)
      {
        auto &dx = icf->dx;
        bfis[dx.vb] += icf->MakeBilinearFormIntegrator ();
      }

    Array<shared_ptr<LinearFormIntegrator>> lfis[4];
    if (lf)
      for (auto icf : lf->icfs)
        {
          auto &dx = icf->dx;
          lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
        }

    vector<shared_ptr<Matrix<SCAL>>> ETmats (ne);
    VVector<SCAL> lfvec (ndof);
    lfvec = 0.0;

    size_t active_elements = 0;
    size_t test_local_ndof = test_fes->GetNDof () / ne;
    Vector<SCAL> sing_val_avg;
    Vector<double> sing_val_max;
    Vector<double> sing_val_min;

    bool has_hidden_dofs = false;
    for (DofId d = 0; d < ndof || !has_hidden_dofs; d++)
      if (HIDDEN_DOF == fes->GetDofCouplingType (d))
        has_hidden_dofs = true;

    ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
      bool definedhere = false;
      for (auto icf : bf->icfs)
        {
          if (icf->dx.vb == VOL)
            if ((!icf->dx.definedonelements)
                || (icf->dx.definedonelements->Test (ei.Nr ())))
              definedhere = true;
        }
      if (!definedhere)
        return; // escape lambda

      Array<DofId> dofs, test_dofs;
      test_fes->GetDofNrs (ei, test_dofs);
      fes->GetDofNrs (ei, dofs);

      FlatMatrix<SCAL> elmat (test_dofs.Size (), dofs.Size (), mlh);

      auto &trafo = ma->GetTrafo (ei, mlh);

      auto &test_fel = test_fes->GetFE (ei, mlh);
      auto &trial_fel = fes->GetFE (ei, mlh);

      elmat = 0.0;
      bool symmetric_so_far = true;
      int bfi_ind = 0;
      while (bfi_ind < bfis[VOL].Size ())
        {
          auto &bfi = bfis[VOL][bfi_ind];
          bfi_ind++;
          if (bfi->DefinedOnElement (ei.Nr ()))
            {
              auto &mapped_trafo
                  = trafo.AddDeformation (bfi->GetDeformation ().get (), mlh);
              try
                {
                  if (mixed_mode)
                    {
                      const auto &mixed_fel
                          = MixedFiniteElement (trial_fel, test_fel);
                      bfi->CalcElementMatrixAdd (mixed_fel, mapped_trafo,
                                                 elmat, symmetric_so_far, mlh);
                    }
                  else
                    {
                      bfi->CalcElementMatrixAdd (test_fel, mapped_trafo, elmat,
                                                 symmetric_so_far, mlh);
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

      if (has_hidden_dofs) // extract visible dofs, if needed
        {
          Array<DofId> vdofs, vtest_dofs;
          fes->GetDofNrs (ei, vdofs, VISIBLE_DOF);
          test_fes->GetDofNrs (ei, vtest_dofs, VISIBLE_DOF);
          FlatMatrix<SCAL> velmat (vtest_dofs.Size (), vdofs.Size (), mlh);
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

      FlatMatrix<SCAL, ColMajor> U (test_dofs.Size (), mlh),
          Vt (dofs.Size (), mlh);
      ngbla::GetSVD<SCAL> (elmat, U, Vt);

      // either tndof is given, or try to determine local size of trefftz
      // space, must be at least dim(trefftz)>=dim(trial)-dim(test)
      int nz = 0;
      if (tndof)
        nz = tndof;
      else
        {
          nz = dofs.Size () - test_dofs.Size ();
          nz = max (nz, 0);
          for (int i = 0; i < min (elmat.Width (), elmat.Height ()); i++)
            if (abs (elmat (i, i)) < eps)
              nz++;
        }
      if (dofs.Size () - test_dofs.Size () > nz)
        throw Exception ("test fes not large enough for given tndof");

      if (getrange)
        ETmats[ei.Nr ()]
            = make_shared<Matrix<SCAL>> (U.Cols (0, dofs.Size () - nz));
      else
        ETmats[ei.Nr ()] = make_shared<Matrix<SCAL>> (
            Trans (Vt.Rows (dofs.Size () - nz, dofs.Size ())));

      if (stats)
        {
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
          int nnz = dofs.Size () - nz;
          auto &test_fel = test_fes->GetFE (ei, mlh);
          auto &trafo = ma->GetTrafo (ei, mlh);
          FlatVector<SCAL> elvec (test_dofs.Size (), mlh),
              elveci (test_dofs.Size (), mlh);
          elvec = 0.0;
          for (auto &lfi : lfis[VOL])
            {
              if (lfi->DefinedOnElement (ei.Nr ()))
                {
                  auto &mapped_trafo = trafo.AddDeformation (
                      lfi->GetDeformation ().get (), mlh);
                  lfi->CalcElementVector (test_fel, mapped_trafo, elveci, mlh);
                  // lfi -> CalcElementVector(fel, mapped_trafo, elveci, mlh);
                  elvec += elveci;
                }
            }
          Matrix<SCAL> Ut = Trans (U).Rows (0, nnz);
          Matrix<SCAL> V = Trans (Vt).Cols (0, nnz);
          Matrix<SCAL> SigI (nnz, nnz);

          SigI = static_cast<SCAL> (0.0);
          // SigI = 0.0;
          for (int i = 0; i < nnz; i++)
            SigI (i, i) = 1.0 / elmat (i, i);
          Matrix<SCAL> elinverse = V * SigI * Ut;

          // lfvec.FV () (dofs) = elinverse * elvec;
          FlatVector<SCAL> elsol (dofs.Size (), mlh);
          elsol = elinverse * elvec;
          lfvec.SetIndirect (dofs, elsol);
        }
    });

    if (stats)
      {
        sing_val_avg /= active_elements;
        (*stats)["singavg"] = Vector<SCAL> (sing_val_avg);
        (*stats)["singmax"] = Vector<double> (sing_val_max);
        (*stats)["singmin"] = Vector<double> (sing_val_min);
      }

    return std::make_tuple (ETmats, make_shared<VVector<SCAL>> (lfvec));
  }

  template std::tuple<vector<shared_ptr<Matrix<double>>>,
                      shared_ptr<BaseVector>>
  EmbTrefftz<double> (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
                      shared_ptr<SumOfIntegrals> lf, double eps,
                      shared_ptr<FESpace> test_fes, int tndof, bool getrange,
                      std::map<std::string, Vector<double>> *stats);
  template std::tuple<vector<shared_ptr<Matrix<Complex>>>,
                      shared_ptr<BaseVector>>
  EmbTrefftz<Complex> (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
                       shared_ptr<SumOfIntegrals> lf, double eps,
                       shared_ptr<FESpace> test_fes, int tndof, bool getrange,
                       std::map<std::string, Vector<Complex>> *stats);

  template <class SCAL>
  shared_ptr<BaseMatrix>
  Elmats2Sparse (vector<shared_ptr<Matrix<SCAL>>> ETmats,
                 shared_ptr<FESpace> fes)
  {
    auto ma = fes->GetMeshAccess ();
    size_t ne = ma->GetNE (VOL);
    size_t ndof = fes->GetNDof ();
    size_t dim = fes->GetDimension ();

    int hidden_dofs = 0;
    for (DofId d : Range (ndof))
      if (HIDDEN_DOF == fes->GetDofCouplingType (d))
        hidden_dofs++;

    Table<int> table, table2;
    TableCreator<int> creator (ne + hidden_dofs);
    TableCreator<int> creator2 (ne + hidden_dofs);
    int prevdofs = 0;

    for (; !creator.Done (); creator++, creator2++)
      {
        prevdofs = 0;
        for (auto ei : ma->Elements (VOL))
          {
            if (!ETmats[ei.Nr ()])
              continue;

            int nz = (ETmats[ei.Nr ()])->Width ();
            Array<DofId> dnums;
            fes->GetDofNrs (ei, dnums, VISIBLE_DOF);
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
                  creator2.Add (ei.Nr (), prevdofs++);
              }
          }

        for (int d = 0, hcnt = 0; d < ndof; d++)
          if (HIDDEN_DOF == fes->GetDofCouplingType (d))
            {
              creator.Add (ne + hcnt, d);
              creator2.Add (ne + hcnt++, prevdofs++);
            }
      }
    table = creator.MoveTable ();
    table2 = creator2.MoveTable ();

    SparseMatrix<SCAL> PP (fes->GetNDof (), prevdofs, table, table2, false);
    auto P = make_shared<SparseMatrix<SCAL>> (PP);
    P->SetZero ();
    for (auto ei : ma->Elements (VOL))
      if (ETmats[ei.Nr ()])
        P->AddElementMatrix (table[ei.Nr ()], table2[ei.Nr ()],
                             *(ETmats[ei.Nr ()]));

    SCAL one = 1;
    FlatMatrix<SCAL> I (1, 1, &one);
    for (int hd = 0; hd < hidden_dofs; hd++)
      P->AddElementMatrix (table[ne + hd], table2[ne + hd], I);

    return P;
  }

  ////////////////////////// EmbTrefftzFESpace ///////////////////////////

  template <typename T, typename shrdT>
  shared_ptr<BaseVector>
  EmbTrefftzFESpace<T, shrdT>::SetOp (shared_ptr<SumOfIntegrals> bf,
                                      shared_ptr<SumOfIntegrals> lf,
                                      double eps, shared_ptr<FESpace> test_fes,
                                      int tndof)
  {
    static Timer timer ("EmbTrefftz: SetOp");

    // shared_ptr<FESpace> use_as_l2 =
    // dynamic_pointer_cast<FESpace>(const_cast<EmbTrefftzFESpace*>(this)->shared_from_this());

    auto embtr = EmbTrefftz<double> (bf, fes, lf, eps, test_fes, tndof, false,
                                     nullptr);
    this->ETmats = std::get<0> (embtr);

    T::Update ();
    int ndof = fes->GetNDof ();
    all2comp.SetSize (ndof);
    all2comp = 0;

    for (auto ei : this->ma->Elements (VOL))
      {
        int nz = (ETmats[ei.Nr ()])->Width ();
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
        this->ctofdof[all2comp[i]] = T::GetDofCouplingType (i);

    T::FinalizeUpdate ();
    return std::get<1> (embtr);
  }

  template <typename T, typename shrdT>
  void EmbTrefftzFESpace<T, shrdT>::GetDofNrs (ElementId ei,
                                               Array<int> &dnums) const
  {
    T::GetDofNrs (ei, dnums);
    if (all2comp.Size () == fes->GetNDof ())
      for (DofId &d : dnums)
        if (IsRegularDof (d))
          d = all2comp[d];
  }

  template <typename T, typename shrdT>
  void EmbTrefftzFESpace<T, shrdT>::VTransformMR (ElementId ei,
                                                  SliceMatrix<double> mat,
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
        mat.Rows (0, nz) = Trans (*(ETmats[ei.Nr ()])) * temp_mat;
      }
  }

  template <typename T, typename shrdT>
  void EmbTrefftzFESpace<T, shrdT>::VTransformVR (ElementId ei,
                                                  SliceVector<double> vec,
                                                  TRANSFORM_TYPE type) const
  {
    static Timer timer ("EmbTrefftz: VTransform");
    RegionTimer reg (timer);

    size_t nz = (ETmats[ei.Nr ()])->Width ();

    if (type == TRANSFORM_RHS)
      {
        Vector<double> new_vec (vec.Size ());
        new_vec = Trans (*(ETmats[ei.Nr ()])) * vec;
        vec = new_vec;
      }
    else if (type == TRANSFORM_SOL)
      {
        Vector<double> new_vec (vec.Size ());
        new_vec = (*(ETmats[ei.Nr ()])) * vec;
        vec = new_vec;
      }
  }

  template class EmbTrefftzFESpace<L2HighOrderFESpace,
                                   shared_ptr<L2HighOrderFESpace>>;
  static RegisterFESpace<
      EmbTrefftzFESpace<L2HighOrderFESpace, shared_ptr<L2HighOrderFESpace>>>
      initembt ("L2EmbTrefftzFESpace");

  // template class EmbTrefftzFESpace<VectorL2FESpace,
  // shared_ptr<VectorL2FESpace>>;
  // static RegisterFESpace<
  // EmbTrefftzFESpace<VectorL2FESpace, shared_ptr<VectorL2FESpace>>>
  // initembt2 ("VL2EmbTrefftzFESpace");

  template class EmbTrefftzFESpace<MonomialFESpace,
                                   shared_ptr<MonomialFESpace>>;
  static RegisterFESpace<
      EmbTrefftzFESpace<MonomialFESpace, shared_ptr<MonomialFESpace>>>
      initembt3 ("MonomialEmbTrefftzFESpace");
}

////////////////////////// python interface ///////////////////////////

#ifdef NGS_PYTHON
template <typename T, typename shrdT>
void ExportETSpace (py::module m, string label)
{
  auto pyspace
      = ngcomp::ExportFESpace<ngcomp::EmbTrefftzFESpace<T, shrdT>> (m, label);

  pyspace.def (py::init ([pyspace] (shrdT fes) {
                 py::list info;
                 auto ma = fes->GetMeshAccess ();
                 info.append (ma);
                 auto nfes
                     = make_shared<ngcomp::EmbTrefftzFESpace<T, shrdT>> (fes);
                 nfes->Update ();
                 nfes->FinalizeUpdate ();
                 connect_auto_update (nfes.get ());
                 return nfes;
               }),
               py::arg ("fes"));

  pyspace.def ("SetOp", &ngcomp::EmbTrefftzFESpace<T, shrdT>::SetOp,
               py::arg ("bf"), py::arg ("lf") = nullptr, py::arg ("eps") = 0,
               py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0);
}

void ExportEmbTrefftz (py::module m)
{
  ExportETSpace<ngcomp::L2HighOrderFESpace,
                shared_ptr<ngcomp::L2HighOrderFESpace>> (
      m, "L2EmbTrefftzFESpace");
  // ExportETSpace<ngcomp::VectorL2FESpace,
  // shared_ptr<ngcomp::VectorL2FESpace>> ( m, "VL2EmbTrefftzFESpace");
  ExportETSpace<ngcomp::MonomialFESpace, shared_ptr<ngcomp::MonomialFESpace>> (
      m, "MonomialEmbTrefftzFESpace");

  m.def (
      "EmbeddedTrefftzFES",
      [] (shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::FESpace> {
        shared_ptr<ngcomp::FESpace> nfes;
        if (dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes))
          nfes = make_shared<ngcomp::EmbTrefftzFESpace<
              ngcomp::L2HighOrderFESpace,
              shared_ptr<ngcomp::L2HighOrderFESpace>>> (
              dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes));
        // else if (dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes))
        // nfes = make_shared<ngcomp::EmbTrefftzFESpace<
        // ngcomp::VectorL2FESpace, shared_ptr<ngcomp::VectorL2FESpace>>> (
        // dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes));
        else if (dynamic_pointer_cast<ngcomp::MonomialFESpace> (fes))
          nfes = make_shared<ngcomp::EmbTrefftzFESpace<
              ngcomp::MonomialFESpace, shared_ptr<ngcomp::MonomialFESpace>>> (
              dynamic_pointer_cast<ngcomp::MonomialFESpace> (fes));
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

  m.def (
      "TrefftzEmbedding",
      [] (shared_ptr<ngfem::SumOfIntegrals> bf,
          shared_ptr<ngcomp::FESpace> fes,
          shared_ptr<ngfem::SumOfIntegrals> lf, double eps,
          shared_ptr<ngcomp::FESpace> test_fes, int tndof, bool getrange,
          optional<py::dict> stats_dict)
          -> std::tuple<shared_ptr<ngcomp::BaseMatrix>,
                        shared_ptr<ngcomp::BaseVector>> {
        shared_ptr<py::dict> pystats = nullptr;
        if (stats_dict)
          pystats = make_shared<py::dict> (*stats_dict);

        if (fes->IsComplex ())
          {
            std::map<std::string, ngcomp::Vector<Complex>> *stats = nullptr;
            if (pystats)
              stats = new std::map<std::string, ngcomp::Vector<Complex>>;
            auto P = ngcomp::EmbTrefftz<Complex> (bf, fes, lf, eps, test_fes,
                                                  tndof, getrange, stats);
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
            auto P = ngcomp::EmbTrefftz<double> (bf, fes, lf, eps, test_fes,
                                                 tndof, getrange, stats);
            if (pystats)
              for (auto const &x : *stats)
                (*pystats)[py::cast (x.first)] = py::cast (x.second);
            return std::make_tuple (
                ngcomp::Elmats2Sparse<double> (std::get<0> (P), fes),
                std::get<1> (P));
          }
      },
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

                :return: [Trefftz embeddint, particular solution]
            )mydelimiter",
      py::arg ("bf"), py::arg ("fes"), py::arg ("lf"), py::arg ("eps") = 0,
      py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0,
      py::arg ("getrange") = false, py::arg ("stats_dict") = nullopt);

  m.def (
      "TrefftzEmbedding",
      [] (shared_ptr<ngfem::SumOfIntegrals> bf,
          shared_ptr<ngcomp::FESpace> fes, double eps,
          shared_ptr<ngcomp::FESpace> test_fes, int tndof, bool getrange,
          optional<py::dict> stats_dict) -> shared_ptr<ngcomp::BaseMatrix> {
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
      },
      R"mydelimiter(
                Used without the parameter lf as input the function only returns the Trefftz embedding.

                :return: Trefftz embedding
            )mydelimiter",
      py::arg ("bf"), py::arg ("fes"), py::arg ("eps") = 0,
      py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0,
      py::arg ("getrange") = false, py::arg ("stats_dict") = py::none ());
}
#endif // NGS_PYTHON
