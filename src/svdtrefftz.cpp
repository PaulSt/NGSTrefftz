#include "svdtrefftz.hpp"

namespace ngcomp
{

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

        dgesdd_ ( &jobz, &m, &n, A.Data(), &lda,
                 S.Data(),
                 U.Data(), &ldu, V.Data(), &ldv,
                 work.Data(), &lwork, iwork.Data(),
                 &info);
        if(info!=0)
            cout << "info = " << info << endl;
        A = 0.0;
        A.Diag(0) = S;
    }

    void LapackSVD (SliceMatrix<double> A,
                    SliceMatrix<double, ColMajor> U,
                    SliceMatrix<double, ColMajor> V)
    {
        Matrix<double,ColMajor> AA = A;
        LapackSVD(AA,U,V);
        A = 0.0;
        A.Diag(0)=AA.Diag();
    }

    std::tuple<shared_ptr<BaseMatrix>,shared_ptr<BaseVector>> SVDTrefftz (shared_ptr<SumOfIntegrals> bf,
                                       shared_ptr<FESpace> fes, double eps,
                                       shared_ptr<SumOfIntegrals> lf
                                       )
    {
        static Timer svdtt("svdtrefftz");
        svdtt.Start();
        LocalHeap lh(1000 * 1000 * 1000);

        auto ma = fes->GetMeshAccess();

        Array<shared_ptr<BilinearFormIntegrator>> bfis[4];  // VOL, BND, ...

        for (auto icf : bf->icfs)
        {
            auto & dx = icf->dx;
            bfis[dx.vb] += make_shared<SymbolicBilinearFormIntegrator> (icf->cf, dx.vb, dx.element_vb);
        }
        Array<shared_ptr<LinearFormIntegrator>> lfis[4];  // VOL, BND, ...
        if(lf)
            for (auto icf : lf->icfs)
            {
                auto & dx = icf->dx;
                lfis[dx.vb] += make_shared<SymbolicLinearFormIntegrator> (icf->cf, dx.vb, dx.element_vb);
            }

        shared_ptr<SparseMatrix<double>> P;
        //Vector<> lfvec(fes->GetNDof());
        VVector<double> lfvec(fes->GetNDof());

        std::once_flag init_flag;
        Table<int> table,table2;

        ma->IterateElements(VOL,lh,[&](auto ei, LocalHeap & mlh)
        {
            HeapReset hr(mlh);
            Array<DofId> dofs;
            fes->GetDofNrs(ei, dofs);
            auto & trafo = ma->GetTrafo(ei, mlh);
            auto & fel = fes->GetFE(ei, mlh);
            FlatMatrix<> elmat(dofs.Size(), mlh), elmati(dofs.Size(), mlh);
            elmat = 0.0;
            for (auto & bfi : bfis[VOL])
            {
                bfi -> CalcElementMatrix(fel, trafo, elmati, mlh);
                elmat += elmati;
            }
            FlatMatrix<double,ColMajor> U(dofs.Size(),mlh), Vt(dofs.Size(),mlh);
            Matrix<> elmato = elmat;
            //FlatMatrix<double,ColMajor> UU(dofs.Size(),mlh), VVt(dofs.Size(),mlh);
            //LapackSVD(elmat,VVt,UU);
            //LapackSVD(elmat,VVt,UU);
            //FlatMatrix<double> U(dofs.Size(),UU.Data()), Vt(dofs.Size(),VVt.Data());
            LapackSVD(elmat,U,Vt);
            //CalcSVD(elmat,U,Vt);
            int nz = 0;
            for(auto sv : elmat.Diag()) if(sv<eps) nz++;


            std::call_once(init_flag, [&](){
                TableCreator<int> creator(ma->GetNE(VOL));
                TableCreator<int> creator2(ma->GetNE(VOL));
                for ( ; !creator.Done(); creator++,creator2++)
                    for(auto ei : ma->Elements(VOL))
                    {
                        Array<DofId> dnums;
                        fes->GetDofNrs (ei, dnums);
                        for (DofId d : dnums)
                            creator.Add (ei.Nr(), d);
                        int prevdofs = (ei.Nr())*nz;
                        for (int d=0;d<nz;d++)
                            creator2.Add (ei.Nr(), d+prevdofs);
                    }
                table = creator.MoveTable();
                table2 = creator2.MoveTable();

                P = make_shared<SparseMatrix<double>>(*(new SparseMatrix<double>(fes->GetNDof(), nz*ma->GetNE(VOL), table,table2,false)));
                P->SetZero();

            });

            //Matrix<> PP = U.Cols(dofs.Size()-nz,dofs.Size());
            Matrix<> PP = Trans(Vt.Rows(dofs.Size()-nz,dofs.Size()));
            P->AddElementMatrix(table[ei.Nr()],table2[ei.Nr()], PP);


            if(lf)
            {
                int nnz = dofs.Size() - nz;
                FlatVector<> elvec(dofs.Size(), mlh), elveci(dofs.Size(), mlh);
                elvec = 0.0;
                for (auto & lfi : lfis[VOL])
                {
                    lfi -> CalcElementVector(fel, trafo, elveci, mlh);
                    elvec += elveci;
                }
                //FlatMatrix<double> Ut(nnz,dofs.Size(),U.Data());
                Matrix<double> Ut = Trans(U).Rows(0,nnz);
                Matrix<double> V = Trans(Vt).Cols(0,nnz);
                //cout << "Vt and V" << endl;
                //cout << Vt << endl;
                //cout << U << endl;
                Matrix<> SigI(nnz,nnz);
                SigI = 0;
                for(int i=0;i<nnz;i++) SigI(i,i)=1.0/elmat(i,i);
                //cout << "elmat, sigi"<<nnz << endl;
                //cout << elmat << endl;
                //cout << SigI << endl;
                //cout << "shapes U,V,Sigi" << endl;
                //cout << Ut.Height() << " " << Ut.Width() << endl;
                //cout << V.Height() << " " << V.Width() << endl;
                //cout << SigI.Height() << " " << SigI.Width() << endl;
                //Matrix<double> calco = U*elmat*VVt;
                //cout << Trans(U*elmat*Vt)<< endl;
                //cout << calco<< endl;
                //cout << "inverse..." << endl;
                //cout << V*SigI*Ut<< endl;
                //cout << "hallo"<<endl;
                //cout << SigI << endl;
                //cout << V << endl;
                //cout << Ut << endl;
                //cout << "hallo"<<endl;
                //cout << Matrix(V*Matrix(SigI*Ut))<< endl;
                Matrix<> elinverse = V*Matrix(SigI*Ut);
            //cout << elinverse*elvec << endl;

                lfvec.FV()(table[ei.Nr()])=elinverse*elvec;

            //cout << (lfvec.FV())(table[ei.Nr()])<<endl;
                //cout << elmato*elinverse<< endl;
                //cout << elmato*elinverse*elmato<< endl;
                //cout << elmato;





            }
        });


        svdtt.Stop();
        //shared_ptr<BaseMatrix> re = make_shared<SparseMatrix<double>>(P);
        return std::make_tuple(P,make_shared<VVector<double>>(lfvec));
    }
}

#ifdef NGS_PYTHON
//#include <python_ngstd.hpp>
//#include <comp.hpp>
//#include <fem.hpp>
//using namespace ngfem;
void ExportSVDTrefftz(py::module m)
{
    m.def("SVDTrefftz", [] (shared_ptr<ngfem::SumOfIntegrals> bf,
                            shared_ptr<ngcomp::FESpace> fes,
                            double eps,
                            shared_ptr<ngfem::SumOfIntegrals> lf
                            )
          {
              return ngcomp::SVDTrefftz(bf,fes,eps,lf);
          },
          py::arg("bf"), py::arg("fes"), py::arg("eps"), py::arg("lf")=nullptr);

    m.def("SVDTrefftz", [] (shared_ptr<ngfem::SumOfIntegrals> bf,
                            shared_ptr<ngcomp::FESpace> fes,
                            double eps
                            ) -> shared_ptr<ngcomp::BaseMatrix>
          {
              return std::get<0>(ngcomp::SVDTrefftz(bf,fes,eps,nullptr));
          },
          py::arg("bf"), py::arg("fes"), py::arg("eps"));

}
#endif // NGS_PYTHON
