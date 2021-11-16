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
        Array<double> work(n*m+100);
        ngbla::integer info;
        char jobu = 'A', jobv = 'A';
        ngbla::integer lda = A.Dist(), ldu = U.Dist(), ldv = V.Dist();
        ngbla::integer lwork = work.Size();

        dgesvd_ ( &jobu, &jobv, &m, &n, A.Data(), &lda,
                 S.Data(),
                 U.Data(), &ldu, V.Data(), &ldv,
                 work.Data(), &lwork,
                 &info);
    if(info!=0)
     cout << "info = " << info << endl;
     //if (n <= 100)
     //cout << "S = " << S << endl;
        A.Diag(0) = S;
    }

    shared_ptr<BaseMatrix> SVDTrefftz (shared_ptr<SumOfIntegrals> bf,
                                       shared_ptr<FESpace> fes, double eps)
    {
        LocalHeap lh(1000 * 1000 * 1000);

        auto ma = fes->GetMeshAccess();

        Array<shared_ptr<BilinearFormIntegrator>> bfis[4];  // VOL, BND, ...

        for (auto icf : bf->icfs)
        {
            auto & dx = icf->dx;
            bfis[dx.vb] += make_shared<SymbolicBilinearFormIntegrator> (icf->cf, dx.vb, dx.element_vb);
        }

        Array<Matrix<double> > elmats;
        Array<Array<int> > newdofs;
        int newdofscounter = 0;
        int max_elsperrow = 0;
        int Pwidth = 0;

        ma->IterateElements(VOL,lh,[&](auto ei, LocalHeap & mlh)
        {
            HeapReset hr(lh);
            Array<DofId> dofs;
            fes->GetDofNrs(ei, dofs);
            auto & trafo = ma->GetTrafo(ei, lh);
            auto & fel = fes->GetFE(ei, lh);

            FlatMatrix<> elmat(dofs.Size(), lh);
            FlatMatrix<> elmati(dofs.Size(), lh);
            elmat = 0.0;
            for (auto & bfi : bfis[VOL])
            {
                bfi -> CalcElementMatrix(fel, trafo, elmati, lh);
                elmat += elmati;
            }
            Matrix<double,ColMajor> U(dofs.Size(),dofs.Size()), V(dofs.Size(),dofs.Size());
            //CalcSVD(elmat,U,V);
            LapackSVD(elmat,U,V);
            Array<int> newdof;
            for(auto sv : elmat.Diag()) if(sv<eps) newdof.Append(newdofscounter++);
            Pwidth += newdof.Size();
            newdofs.Append(newdof);
            elmats.Append(U.Cols(dofs.Size()-newdof.Size(),dofs.Size()));
        });

        int nnewndof = newdofs.Last().Last();


        TableCreator<int> creator(ma->GetNE(VOL));
        TableCreator<int> creator2(ma->GetNE(VOL));
        for ( ; !creator.Done(); creator++, creator2++)
            ma->IterateElements(VOL,lh,[&](auto ei, LocalHeap & mlh)
                                {
                                    Array<DofId> dnums;
                                    fes->GetDofNrs (ei, dnums);
                                    for (DofId d : dnums)
                                        if (IsRegularDof(d)) creator.Add (ei.Nr(), d);
                                    for (int d : newdofs[ei.Nr()])
                                        if (IsRegularDof(d)) creator2.Add (ei.Nr(), d);
                                });
        auto table = creator.MoveTable();
        auto table2 = creator2.MoveTable();

        SparseMatrix<double> P(fes->GetNDof(), Pwidth, table,table2,false);

        for(auto ei : ma->Elements())
        {
            Array<DofId> dofs;
            fes->GetDofNrs(ei, dofs);
            P.AddElementMatrix(dofs,newdofs[ei.Nr()], elmats[ei.Nr()]);
        }
        return make_shared<SparseMatrix<double>>(P);
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
              return SVDTrefftz(bf,fes,eps);
          },
          py::arg("bf"), py::arg("fes"), py::arg("eps"));
}
#endif // NGS_PYTHON
