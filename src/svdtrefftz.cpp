#include "svdtrefftz.hpp"

namespace ngcomp
{

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

        //Array<double> allvalues(fes->GetNDof()*fes->GetNDof());
        Array<Matrix<double> > Pmats(ma->GetNE(VOL));
        Array<int> nzs(ma->GetNE(VOL));

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
            CalcSVD(elmat,U,V);
            int nz = 0;
            for(auto sv : elmat.Diag()) if(sv<eps) nz++;
            nzs[ei.Nr()]=nz;
            //Pmats[ei.Nr()].AssignMemory (dofs.Size(),dofs.Size(), allvalues.Addr(dofs.Size()*dofs.Size()*ei.Nr()));
            Pmats[ei.Nr()] = (U.Cols(dofs.Size()-nz,dofs.Size()));
        }
        );

        TableCreator<int> creator(ma->GetNE(VOL));
        TableCreator<int> creator2(ma->GetNE(VOL));
        for ( ; !creator.Done(); creator++,creator2++)
            ma->IterateElements(VOL,[&](auto ei)
                                {
                                    Array<DofId> dnums;
                                    fes->GetDofNrs (ei, dnums);
                                    for (DofId d : dnums)
                                        //if (IsRegularDof(d))
                                        creator.Add (ei.Nr(), d);
                                    int prevdofs = 0;
                                    for (int i=0;i<ei.Nr();i++) prevdofs += nzs[ei.Nr()];
                                    for (int d=0;d<nzs[ei.Nr()];d++)
                                        creator2.Add (ei.Nr(), d+prevdofs);
                                }
        );
        auto table = creator.MoveTable();
        auto table2 = creator2.MoveTable();

        int Pwidth = 0;
        for(auto ei : ma->Elements(VOL)) Pwidth += nzs[ei.Nr()];
        if(Pwidth != nzs[0]*(ma->GetNE(VOL))) throw std::logic_error("WARNING: inconsistent nzs, try larger eps?");

        SparseMatrix<double> P(fes->GetNDof(), Pwidth, table,table2,false);
        P.SetZero();

        ma->IterateElements(VOL,[&](auto ei)
        {
            P.AddElementMatrix(table[ei.Nr()],table2[ei.Nr()], Pmats[ei.Nr()]);
        }
        );
        svdtt.Stop();
        shared_ptr<BaseMatrix> re = make_shared<SparseMatrix<double>>(P);
        return re;
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
