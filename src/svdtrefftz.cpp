#include "svdtrefftz.hpp"

namespace ngcomp
{

    void LapackSVD (SliceMatrix<> A,
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

    shared_ptr<BaseMatrix> SVDTrefftz (shared_ptr<SumOfIntegrals> bf,
                                       shared_ptr<FESpace> fes, double eps)
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

        shared_ptr<SparseMatrix<double>> P;
        Table<int> table,table2;

        std::once_flag init_flag;

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
            FlatMatrix<double,ColMajor> U(dofs.Size(),mlh), V(dofs.Size(),mlh);
            LapackSVD(elmat,U,V);
            //CalcSVD(elmat,U,V);
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

            Matrix<> PP = U.Cols(dofs.Size()-nz,dofs.Size());
            P->AddElementMatrix(table[ei.Nr()],table2[ei.Nr()], PP);
        });


        svdtt.Stop();
        return P;
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
                            double eps) -> shared_ptr<ngcomp::BaseMatrix>
          {
              return ngcomp::SVDTrefftz(bf,fes,eps);
          },
          py::arg("bf"), py::arg("fes"), py::arg("eps"));

}
#endif // NGS_PYTHON
