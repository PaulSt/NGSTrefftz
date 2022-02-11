#include "embtrefftz.hpp"
#include <bla.hpp>

namespace ngbla{

#ifdef LAPACK
    void LapackSVD (SliceMatrix<double, ColMajor> A,
                    SliceMatrix<double, ColMajor> U,
                    SliceMatrix<double, ColMajor> V)
    {
        static Timer t("LapackSVD"); RegionTimer reg(t);
        ngbla::integer m = A.Width(), n = A.Height();
        Vector<> S(min(n,m));
        Array<double> work(4*m*m+6*m+m+100);
        Array<int> iwork(max(n,m)*9);
        ngbla::integer info;
        char jobz = 'A';
        ngbla::integer lda = A.Dist(), ldu = U.Dist(), ldv = V.Dist();
        ngbla::integer lwork = work.Size();

        //if(std::is_same<SCAL,double>::value)
        dgesdd_ ( &jobz, &m, &n, A.Data(), &lda,
                 S.Data(),
                 U.Data(), &ldu, V.Data(), &ldv,
                 work.Data(), &lwork, iwork.Data(),
                 &info);
        if(info!=0)
            throw Exception("something went wrong in the svd " + std::to_string(info));
        A = 0.0;
        A.Diag(0) = S;
    }

    void LapackSVD (SliceMatrix<Complex, ColMajor> A,
                    SliceMatrix<Complex, ColMajor> U,
                    SliceMatrix<Complex, ColMajor> V)
    {
        static Timer t("LapackSVD"); RegionTimer reg(t);
        ngbla::integer m = A.Width(), n = A.Height();
        Vector<> S(min(n,m));
        Array<Complex> work(4*m*m+6*m+m+100);
        Array<int> iwork(max(n,m)*9);
        Array<double> rwork(5*max(n,m)*min(n,m)+5*min(n,m));
        ngbla::integer info;
        char jobz = 'A';
        ngbla::integer lda = A.Dist(), ldu = U.Dist(), ldv = V.Dist();
        ngbla::integer lwork = work.Size();

        zgesdd_ ( &jobz, &m, &n, A.Data(), &lda,
                 S.Data(),
                 U.Data(), &ldu, V.Data(), &ldv,
                 work.Data(), &lwork, rwork.Data(), iwork.Data(),
                 &info);
        if(info!=0)
            throw Exception("something went wrong in the svd " + std::to_string(info));
        A = 0.0;
        A.Diag(0) = S;
    }
#endif

    template <class SCAL>
    void GetSVD (SliceMatrix<SCAL> A,
                    SliceMatrix<SCAL, ColMajor> U,
                    SliceMatrix<SCAL, ColMajor> V)
    {
        Matrix<SCAL,ColMajor> AA = A;
        //Matrix<SCAL,ColMajor> AA(A.Height(),A.Width());
        //for(int i=0;i<A.Height();i++)
            //for(int j=0;j<A.Width();j++)
                //AA(i,j)= A(i,j);
#ifdef LAPACK
        LapackSVD(AA,U,V);
#else
        CalcSVD(AA,U,V);
#endif
        A = 0.0;
        A.Diag(0)=AA.Diag();
    }

    template
    void GetSVD<double>
        (SliceMatrix<double> A, SliceMatrix<double, ColMajor> U, SliceMatrix<double, ColMajor> V);

    template
    void GetSVD<Complex>
        (SliceMatrix<Complex> A, SliceMatrix<Complex, ColMajor> U, SliceMatrix<Complex, ColMajor> V);
}


namespace ngcomp
{

    template <class SCAL>
    std::tuple<shared_ptr<BaseMatrix>,shared_ptr<BaseVector>> EmbTrefftz (shared_ptr<SumOfIntegrals> bf,
                                       shared_ptr<FESpace> fes, double eps,
                                       shared_ptr<SumOfIntegrals> lf
                                       )
    {
        static Timer svdtt("svdtrefftz"); RegionTimer reg(svdtt);
        LocalHeap lh(1000 * 1000 * 1000);

        auto ma = fes->GetMeshAccess();

        Array<shared_ptr<BilinearFormIntegrator>> bfis[4];  // VOL, BND, ...
        for (auto icf : bf->icfs)
        {
            auto & dx = icf->dx;
            bfis[dx.vb] += icf->MakeBilinearFormIntegrator();
        }

        Array<shared_ptr<LinearFormIntegrator>> lfis[4];
        if(lf)
            for (auto icf : lf->icfs)
            {
                auto & dx = icf->dx;
                lfis[dx.vb] += icf->MakeLinearFormIntegrator();
            }

        shared_ptr<SparseMatrix<SCAL>> P;
        VVector<SCAL> lfvec(fes->GetNDof());

        std::once_flag init_flag;
        Table<int> table,table2;

        ma->IterateElements(VOL,lh,[&](auto ei, LocalHeap & mlh)
        {
            HeapReset hr(mlh);
            Array<DofId> dofs;
            fes->GetDofNrs(ei, dofs);

            bool definedhere = false;
            for (auto icf : bf->icfs)
            {
                if (icf->dx.vb == VOL)
                    if ( (!icf->dx.definedonelements) || (icf->dx.definedonelements->Test(ei.Nr())))
                        definedhere = true;
            }
            if (!definedhere)
                return; // escape lambda

            auto & trafo = ma->GetTrafo(ei, mlh);
            auto & fel = fes->GetFE(ei, mlh);
            
            FlatMatrix<SCAL> elmat(dofs.Size(), mlh);
            elmat = 0.0;
            bool symmetric_so_far = true;

            int bfi_ind = 0;
            while (bfi_ind < bfis[VOL].Size())
            {
                auto & bfi = bfis[VOL][bfi_ind];
                bfi_ind++;
                if (bfi->DefinedOnElement(ei.Nr()))
                {
                    auto & mapped_trafo = trafo.AddDeformation(bfi->GetDeformation().get(), mlh);
                    try
                    {
                        bfi -> CalcElementMatrixAdd(fel, mapped_trafo, elmat, symmetric_so_far, mlh);
                    }
                    catch (ExceptionNOSIMD e)
                    {
                        elmat = 0.0;
                        cout << IM(6) << "ExceptionNOSIMD " << e.What() << endl
                        << "switching to scalar evaluation" << endl;
                        bfi -> SetSimdEvaluate (false);
                        bfi_ind = 0;
                    }
                }
            }
            FlatMatrix<SCAL,ColMajor> U(dofs.Size(),mlh), Vt(dofs.Size(),mlh);
            ngbla::GetSVD<SCAL>(elmat,U,Vt);

            int nz = 0;
            // assumption here: all (active) elements have the same number of (weak) Trefftz fcts.
            for(auto sv : elmat.Diag()) if(abs(sv)<eps) nz++;

            std::call_once(init_flag, [&](){
                TableCreator<int> creator(ma->GetNE(VOL));
                TableCreator<int> creator2(ma->GetNE(VOL));
                int prevdofs = 0;
                for ( ; !creator.Done(); creator++,creator2++)
                {
                    prevdofs = 0;
                    for(auto ei : ma->Elements(VOL))
                    {
                        Array<DofId> dnums;
                        fes->GetDofNrs (ei, dnums);
                        bool hasregdof = false;
                        for (DofId d : dnums)
                            if (IsRegularDof(d))
                            {
                                creator.Add (ei.Nr(), d);
                                hasregdof = true;
                            }
                        // assumption here: Either all or no dof is regular
                        if (hasregdof)
                        {
                            for (int d=0;d<nz;d++)
                                creator2.Add (ei.Nr(), d+prevdofs);
                            prevdofs += nz;
                        }
                    }
                }
                table = creator.MoveTable();
                table2 = creator2.MoveTable();

                SparseMatrix<SCAL> PP(fes->GetNDof(), prevdofs, table,table2,false);
                P = make_shared<SparseMatrix<SCAL>>(PP);
                P->SetZero();
            });

            Matrix<SCAL> PP = Trans(Vt.Rows(dofs.Size()-nz,dofs.Size()));
            P->AddElementMatrix(table[ei.Nr()],table2[ei.Nr()], PP);


            if(lf)
            {
                int nnz = dofs.Size() - nz;
                FlatVector<SCAL> elvec(dofs.Size(), mlh), elveci(dofs.Size(), mlh);
                elvec = 0.0;
                for (auto & lfi : lfis[VOL])
                {
                    if (lfi->DefinedOnElement(ei.Nr()))
                    {
                        auto & mapped_trafo = trafo.AddDeformation(lfi->GetDeformation().get(), mlh);
                        lfi -> CalcElementVector(fel, mapped_trafo, elveci, mlh);
                        elvec += elveci;
                    }
                }
                Matrix<SCAL> Ut = Trans(U).Rows(0,nnz);
                Matrix<SCAL> V = Trans(Vt).Cols(0,nnz);
                Matrix<SCAL> SigI(nnz,nnz);
                SigI = static_cast<SCAL>(0.0);
                //SigI = 0.0;
                for(int i=0;i<nnz;i++) SigI(i,i)=1.0/elmat(i,i);
                Matrix<SCAL> elinverse = V*SigI*Ut;
                lfvec.FV()(table[ei.Nr()])=elinverse*elvec;
            }
        });

        return std::make_tuple(P,make_shared<VVector<SCAL>>(lfvec));
    }

  template
      std::tuple<shared_ptr<BaseMatrix>,shared_ptr<BaseVector>> EmbTrefftz<double>
          (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes, double eps, shared_ptr<SumOfIntegrals> lf);
  template
      std::tuple<shared_ptr<BaseMatrix>,shared_ptr<BaseVector>> EmbTrefftz<Complex>
          (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes, double eps, shared_ptr<SumOfIntegrals> lf);

}

#ifdef NGS_PYTHON
//#include <python_ngstd.hpp>
//#include <comp.hpp>
//#include <fem.hpp>
//using namespace ngfem;
void ExportEmbTrefftz(py::module m)
{
    m.def("TrefftzEmbedding", [] (shared_ptr<ngfem::SumOfIntegrals> bf,
                            shared_ptr<ngcomp::FESpace> fes,
                            double eps,
                            shared_ptr<ngfem::SumOfIntegrals> lf
                            )
          {
              if(fes->IsComplex())
                  return ngcomp::EmbTrefftz<Complex>(bf,fes,eps,nullptr);

              return ngcomp::EmbTrefftz<double>(bf,fes,eps,lf);
          },
          py::arg("bf"), py::arg("fes"), py::arg("eps"), py::arg("lf")=nullptr);

    m.def("TrefftzEmbedding", [] (shared_ptr<ngfem::SumOfIntegrals> bf,
                            shared_ptr<ngcomp::FESpace> fes,
                            double eps
                            ) -> shared_ptr<ngcomp::BaseMatrix>
          {
              if(fes->IsComplex())
                  return std::get<0>(ngcomp::EmbTrefftz<Complex>(bf,fes,eps,nullptr));

              return std::get<0>(ngcomp::EmbTrefftz<double>(bf,fes,eps,nullptr));
          },
          py::arg("bf"), py::arg("fes"), py::arg("eps"));

}
#endif // NGS_PYTHON
