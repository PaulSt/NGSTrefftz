#include "svdtrefftz.hpp"
//#include "../fem/integratorcf.hpp"
//#include "../fem/h1lofe.hpp"


namespace ngcomp
{
    shared_ptr<BaseMatrix> SVDTrefftz (shared_ptr<SumOfIntegrals> bf,
                     shared_ptr<FESpace> fes, double eps)
    {
        LocalHeap lh(1000 * 1000 * 1000);

        auto ma = fes->GetMeshAccess();
        //cout << *bf;

        //const BitArray & freedofs = *fes->GetFreeDofs();

        Array<shared_ptr<BilinearFormIntegrator>> bfis[4];  // VOL, BND, ...
        //auto gfvec = gf->GetVector().FV<double>();
        //gfvec = 0.0;

        for (auto icf : bf->icfs)
        {
            auto & dx = icf->dx;
            bfis[dx.vb] += make_shared<SymbolicBilinearFormIntegrator> (icf->cf, dx.vb, dx.element_vb);
        }

        //Array<FlatMatrix<double> > elmats(ma->GetNE());
        Array<Matrix<double> > elmats;
        Array<int> widths;
        size_t Tndofs = 0;

        ma->IterateElements(VOL,lh,[&](auto ei, LocalHeap & mlh)
            //for(auto ei : ma->Elements())
            {
                HeapReset hr(lh);
                //ElementId ei(VOL, el);
                //cout << "ei = " << ei << endl;
                Array<DofId> dofs;
                fes->GetDofNrs(ei, dofs);
                //cout << "dofs = " << dofs << endl;
                auto & trafo = ma->GetTrafo(ei, lh);
                auto & fel = fes->GetFE(ei, lh);

                FlatMatrix<> elmat(dofs.Size(), lh);
                FlatMatrix<> elmati(dofs.Size(), lh);
                elmat = 0.0;
                for (auto & bfi : bfis[VOL])
                {
                    bfi -> CalcElementMatrix(fel, trafo, elmati, lh);
                    elmat += elmati;
                    // bfi -> CalcElementMatrixAdd(fel, trafo, elmat, lh);
                }
                //cout << ei << " elmat = " << elmat << endl;
                SliceMatrix<> A(elmat);
                Matrix<double,ColMajor> U(dofs.Size(),dofs.Size()), V(dofs.Size(),dofs.Size());
                CalcSVD(elmat,U,V);
                int nnz=0;
                for(auto sv : elmat.Diag()) if(sv>eps) nnz++;
                //cout << nnz <<  "SVD" << elmat.Diag() << endl << U << endl << V << endl;
                //FlatMatrix<> P = U.Cols(nnz,dofs.Size());
                //cout <<"P"<<endl<< P << endl;
                elmats.Append(U.Cols(nnz,dofs.Size()));
                widths.Append(dofs.Size()-nnz);
                Tndofs += dofs.Size()-nnz;
                //(elmats[ei.Nr()]
            });

        ElementByElementMatrix<double> P(fes->GetNDof(), Tndofs, ma->GetNE(), false);
        for(auto ei : ma->Elements())
            {
                Array<DofId> dofs;
                fes->GetDofNrs(ei, dofs);
                P.AddElementMatrix(ei.Nr(),dofs,dofs.Range(0,widths[ei.Nr()]), elmats[ei.Nr()]);
            }
        return make_shared<ElementByElementMatrix<double>>(P);
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
                            shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::BaseMatrix>
          {
              return SVDTrefftz(bf,fes);
          },
          py::arg("bf"), py::arg("fes"));
}
#endif // NGS_PYTHON
