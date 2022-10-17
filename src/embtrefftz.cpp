#include "embtrefftz.hpp"

namespace ngbla
{

#ifdef LAPACK
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
#ifdef LAPACK
    LapackSVD (AA, U, V);
#else
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

  template void
  GetSVD<Complex> (SliceMatrix<Complex> A, SliceMatrix<Complex, ColMajor> U,
                   SliceMatrix<Complex, ColMajor> V);
}

namespace ngcomp
{

  template <class SCAL>
  std::tuple<shared_ptr<Array<Matrix<SCAL>>>, shared_ptr<BaseVector>>
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

    // shared_ptr<SparseMatrix<SCAL>> P;

    Array<Matrix<SCAL>> ETmats (ne);
    VVector<SCAL> lfvec (ndof);

    // std::once_flag init_flag;
    // Table<int> table,table2;

    Array<Vector<SCAL>> singular_values (ma->GetNE ());
    size_t active_elements = 0;

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
      active_elements += 1;

      Array<DofId> test_dofs;
      test_fes->GetDofNrs (ei, test_dofs);
      Array<DofId> dofs;
      fes->GetDofNrs (ei, dofs);

      FlatMatrix<SCAL> elmat (test_dofs.Size (), dofs.Size (), mlh);
      FlatMatrix<SCAL, ColMajor> U (test_dofs.Size (), mlh),
          Vt (dofs.Size (), mlh);

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
      ngbla::GetSVD<SCAL> (elmat, U, Vt);

      // assumption here: all (active) elements have the same number of (weak)
      // Trefftz fcts.
      int nz = 0;
      if (tndof)
        nz = tndof;
      else
        {
          nz = trial_fel.GetNDof () - test_fel.GetNDof ();
          for (int i = 0; i < min (elmat.Width (), elmat.Height ()); i++)
            if (abs (elmat (i, i)) < eps)
              nz++;
        }

      if (stats)
        singular_values[ei.Nr ()] = Vector<SCAL> (elmat.Diag ());

      // Matrix<SCAL> PP;
      if (getrange)
        ETmats[ei.Nr ()] = U.Cols (0, dofs.Size () - nz);
      else
        ETmats[ei.Nr ()] = Trans (Vt.Rows (dofs.Size () - nz, dofs.Size ()));
      // P->AddElementMatrix(table[ei.Nr()],table2[ei.Nr()], PP);

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

          lfvec.FV () (dofs) = elinverse * elvec;
        }
    });

    if (stats)
      {
        // Singular value avg
        size_t maxdofs = 0;
        for (auto v : singular_values)
          maxdofs = max (maxdofs, v.Size ());
        Vector<SCAL> sing_val_avg (maxdofs);
        sing_val_avg = 0;
        for (auto v : singular_values)
          if (v.Size () > 0)
            sing_val_avg += v;
        sing_val_avg /= active_elements;
        (*stats)["singavg"] = Vector<SCAL> (sing_val_avg);

        // Singular value variance
        Vector<SCAL> sing_val_var (maxdofs);
        sing_val_var = 0;
        for (auto v : singular_values)
          for (int j = 0; j < v.Size (); j++)
            sing_val_var[j] += pow (v[j] - sing_val_avg[j], 2);
        sing_val_var /= active_elements;
        for (auto &v : sing_val_var)
          v = sqrt (v);
        (*stats)["singvar"] = Vector<SCAL> (sing_val_var);
      }

    return std::make_tuple (make_shared<Array<Matrix<SCAL>>> (ETmats),
                            make_shared<VVector<SCAL>> (lfvec));
  }

  template std::tuple<shared_ptr<Array<Matrix<double>>>,
                      shared_ptr<BaseVector>>
  EmbTrefftz<double> (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
                      shared_ptr<SumOfIntegrals> lf, double eps,
                      shared_ptr<FESpace> test_fes, int tndof, bool getrange,
                      std::map<std::string, Vector<double>> *stats);
  template std::tuple<shared_ptr<Array<Matrix<Complex>>>,
                      shared_ptr<BaseVector>>
  EmbTrefftz<Complex> (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
                       shared_ptr<SumOfIntegrals> lf, double eps,
                       shared_ptr<FESpace> test_fes, int tndof, bool getrange,
                       std::map<std::string, Vector<Complex>> *stats);

  template <class SCAL>
  shared_ptr<BaseMatrix> Elmats2Sparse (shared_ptr<Array<Matrix<SCAL>>> ETmats,
                                        shared_ptr<FESpace> fes)
  {
    auto ma = fes->GetMeshAccess ();
    size_t ne = ma->GetNE (VOL);
    size_t ndof = fes->GetNDof ();
    size_t dim = fes->GetDimension ();

    Table<int> table, table2;
    TableCreator<int> creator (ne);
    TableCreator<int> creator2 (ne);
    int prevdofs = 0;
    for (; !creator.Done (); creator++, creator2++)
      {
        prevdofs = 0;
        for (auto ei : ma->Elements (VOL))
          {
            int nz = (*ETmats)[ei.Nr ()].Width ();
            Array<DofId> dnums;
            fes->GetDofNrs (ei, dnums);
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
                  creator2.Add (ei.Nr (), d + prevdofs);
                prevdofs += nz;
              }
          }
      }
    table = creator.MoveTable ();
    table2 = creator2.MoveTable ();

    SparseMatrix<SCAL> PP (fes->GetNDof (), prevdofs, table, table2, false);
    auto P = make_shared<SparseMatrix<SCAL>> (PP);
    P->SetZero ();
    for (auto ei : ma->Elements (VOL))
      P->AddElementMatrix (table[ei.Nr ()], table2[ei.Nr ()],
                           (*ETmats)[ei.Nr ()]);
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
        int nz = (*ETmats)[ei.Nr ()].Width ();
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

  template class EmbTrefftzFESpace<L2HighOrderFESpace,
                                   shared_ptr<L2HighOrderFESpace>>;
  static RegisterFESpace<
      EmbTrefftzFESpace<L2HighOrderFESpace, shared_ptr<L2HighOrderFESpace>>>
      initembt ("L2EmbTrefftzFESpace");
  template class EmbTrefftzFESpace<VectorL2FESpace,
                                   shared_ptr<VectorL2FESpace>>;
  static RegisterFESpace<
      EmbTrefftzFESpace<VectorL2FESpace, shared_ptr<VectorL2FESpace>>>
      initembt2 ("VL2EmbTrefftzFESpace");
}

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
  ExportETSpace<ngcomp::VectorL2FESpace, shared_ptr<ngcomp::VectorL2FESpace>> (
      m, "VL2EmbTrefftzFESpace");

  m.def (
      "EmbeddedTrefftzFES",
      [] (shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::FESpace> {
        shared_ptr<ngcomp::FESpace> nfes;
        if (dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes))
          nfes = make_shared<ngcomp::EmbTrefftzFESpace<
              ngcomp::L2HighOrderFESpace,
              shared_ptr<ngcomp::L2HighOrderFESpace>>> (
              dynamic_pointer_cast<ngcomp::L2HighOrderFESpace> (fes));
        else if (dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes))
          nfes = make_shared<ngcomp::EmbTrefftzFESpace<
              ngcomp::VectorL2FESpace, shared_ptr<ngcomp::VectorL2FESpace>>> (
              dynamic_pointer_cast<ngcomp::VectorL2FESpace> (fes));
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
          py::object stats_dict)
          -> std::tuple<shared_ptr<ngcomp::BaseMatrix>,
                        shared_ptr<ngcomp::BaseVector>> {
        py::extract<py::dict> stats_ (stats_dict);
        shared_ptr<py::dict> pystats = nullptr;
        if (stats_.check ())
          pystats = make_shared<py::dict> (stats_ ());

        if (fes->IsComplex ())
          {
            std::map<std::string, ngcomp::Vector<Complex>> stats;
            auto P = ngcomp::EmbTrefftz<Complex> (bf, fes, lf, eps, test_fes,
                                                  tndof, getrange, &stats);
            if (pystats)
              for (auto const &x : stats)
                (*pystats)[py::cast (x.first)] = py::cast (x.second);
            return std::make_tuple (
                ngcomp::Elmats2Sparse<Complex> (std::get<0> (P), fes),
                std::get<1> (P));
          }
        else
          {
            std::map<std::string, ngcomp::Vector<double>> stats;
            auto P = ngcomp::EmbTrefftz<double> (bf, fes, lf, eps, test_fes,
                                                 tndof, getrange, &stats);
            if (pystats)
              for (auto const &x : stats)
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
                :param tndof: If known, local ndofs of the Trefftz space, also eps and/or test_fes are used to find the dimension, defaults to 0
                :param getrange: If True, extract the range instead of the kernel
                :param stats_dict: Pass a dictionary to fill it with stats on the singular values.

                :return: [Trefftz embeddint, particular solution]
            )mydelimiter",
      py::arg ("bf"), py::arg ("fes"), py::arg ("lf"), py::arg ("eps") = 0,
      py::arg ("test_fes") = nullptr, py::arg ("tndof") = 0,
      py::arg ("getrange") = false, py::arg ("stats_dict") = py::none ());

  m.def (
      "TrefftzEmbedding",
      [] (shared_ptr<ngfem::SumOfIntegrals> bf,
          shared_ptr<ngcomp::FESpace> fes, double eps,
          shared_ptr<ngcomp::FESpace> test_fes, int tndof, bool getrange,
          py::object stats_dict) -> shared_ptr<ngcomp::BaseMatrix> {
        py::extract<py::dict> stats_ (stats_dict);
        shared_ptr<py::dict> pystats = nullptr;
        if (stats_.check ())
          pystats = make_shared<py::dict> (stats_ ());

        if (fes->IsComplex ())
          {
            std::map<std::string, ngcomp::Vector<Complex>> stats;
            auto P = std::get<0> (ngcomp::EmbTrefftz<Complex> (
                bf, fes, nullptr, eps, test_fes, tndof, getrange, &stats));
            if (pystats)
              for (auto const &x : stats)
                (*pystats)[py::cast (x.first)] = py::cast (x.second);
            return ngcomp::Elmats2Sparse<Complex> (P, fes);
          }
        else
          {
            std::map<std::string, ngcomp::Vector<double>> stats;
            auto P = std::get<0> (ngcomp::EmbTrefftz<double> (
                bf, fes, nullptr, eps, test_fes, tndof, getrange, &stats));
            if (pystats)
              for (auto const &x : stats)
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
