#include "svdtrefftz.hpp"
//#include "../fem/integratorcf.hpp"
//#include "../fem/h1lofe.hpp"


namespace ngcomp
{

    class ConstCOO
    {
        public:
            ConstCOO(){;}
            Array<int> rowind;
            Array<int> colind;
            Array<double> data;

            void AppendElmat(const Array<int> & dofs1,const Array<int> & dofs2, FlatMatrix<> elmat)
            {
                for(int i=0;i<dofs1.Size();i++)
                    for(int j=0;j<dofs2.Size();j++)
                    {
                        rowind.Append(dofs1[i]);
                        colind.Append(dofs2[j]);
                        data.Append(elmat(i,j));
                    }
            };
    };

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
        //Array<Matrix<double> > elmats;
        Array<Array<int> > newdofs;
        int newdofscounter = 0;
        //Array<int> widths;
        //size_t Tndofs = 0;
        int max_elsperrow = 0;
        int Pwidth = 0;
        ConstCOO elmats;

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
                //int nnz=0;
                Array<int> newdof;
                for(auto sv : elmat.Diag()) if(sv<eps) newdof.Append(newdofscounter++); 
                Pwidth += newdof.Size();
                //newdofs.Append(newdof);
                //max_elsperrow = max_elsperrow<newdof.Size()?newdof.Size():max_elsperrow;
                //cout << newdof.Size();
                //cout << nnz <<  "SVD" << elmat.Diag() << endl << U << endl << V << endl;
                //FlatMatrix<> P = U.Cols(nnz,dofs.Size());
                //cout <<"P"<<endl<< P << endl;
                //elmats.Append(U.Cols(dofs.Size()-newdof.Size(),dofs.Size()));
                //widths.Append(dofs.Size()-nnz);
                //Tndofs += dofs.Size()-nnz;
                //(elmats[ei.Nr()]
                elmats.AppendElmat(dofs,newdof,elmat);
            });

        //int nnewndof = newdofs.Last().Last();
        //cout << endl<<nnewndof<<endl;
        //ElementByElementMatrix<double> P(fes->GetNDof(), 29, ma->GetNE(), false);
        //cout << fes->GetNDof() << ","<<max_elsperrow << endl;
        //SparseMatrix<double> P(fes->GetNDof(), max_elsperrow);
        //cout << "fisti"<<P.First(0)<<endl;
        //cout << "fisti"<<P.GetFirstArray()<<endl;
        //cout << "colin"<<P.GetColIndices()<<endl;
        //cout << P.GetRowIndices(0) << endl;
        //for(auto ei : ma->Elements())
            //{
                //Array<DofId> dofs;
                //fes->GetDofNrs(ei, dofs);
                ////cout << dofs << " and " << newdofs[ei.Nr()] << endl;
                ////P.AddElementMatrix(ei.Nr(),dofs,newdofs[ei.Nr()], elmats[ei.Nr()]);
                //for(int i : dofs)
                    //for(int j : newdofs[ei.Nr()])
                    //{
                        //P.CreatePosition(i,j);
                        //cout << i<<","<<j<< endl;
                    //}
                //P.AddElementMatrix(dofs,newdofs[ei.Nr()], elmats[ei.Nr()]);
                ////cout << dofs.Range(0,widths[ei.Nr()]) << endl;
                ////cout << newdofs[ei.Nr()];
            //}
        //return make_shared<ElementByElementMatrix<double>>(P);
        
        shared_ptr<SparseMatrixTM<double>> P = SparseMatrixTM<double>::CreateFromCOO(elmats.rowind,elmats.colind,elmats.data,fes->GetNDof(),Pwidth);
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
                            shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<ngcomp::BaseMatrix>
          {
              return SVDTrefftz(bf,fes);
          },
          py::arg("bf"), py::arg("fes"));
}
#endif // NGS_PYTHON
