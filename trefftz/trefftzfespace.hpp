#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP
#include "trefftzwavefe.hpp"


namespace ngcomp
{

    class TrefftzFESpace : public FESpace
    {
        int D;
        int fullD;
        int order;
        size_t ndof;
        int nel;
        int nvert;
        int local_ndof;
        float c=1;
        int useshift=1;
        int usescale=1;
        int useqt=0;
        int basistype;
        shared_ptr<CoefficientFunction> wavespeedcf;
        Matrix<shared_ptr<CoefficientFunction>> wavespeedmatrix;
            Array<Matrix<double>> gamma;

        public:


        TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

        void SetWavespeed(shared_ptr<CoefficientFunction> awavespeedcf) {
            wavespeedcf=awavespeedcf;
            if(useqt)
            {
                cout << "started auto diff.... ";
                shared_ptr<CoefficientFunction> localwavespeedcf = make_shared<ConstantCoefficientFunction>(1)/(wavespeedcf*wavespeedcf);
                shared_ptr<CoefficientFunction> localwavespeedcfx = make_shared<ConstantCoefficientFunction>(1)/(wavespeedcf*wavespeedcf);
                //wavespeedmatrix.SetSize(this->order-1);
                //for(int ny=0;ny<=(this->order-2)*(D==2);ny++)
                //{
                //for(int nx=0;nx<this->order-1;nx++)
                //{
                //wavespeedmatrix(nx,ny) = localwavespeedcfx;
                //localwavespeedcfx = localwavespeedcfx->Diff(MakeCoordinateCoefficientFunction(0).get(), make_shared<ConstantCoefficientFunction>(1) );
                //}
                //localwavespeedcf = localwavespeedcf->Diff(MakeCoordinateCoefficientFunction(1).get(), make_shared<ConstantCoefficientFunction>(1) );
                //localwavespeedcfx = localwavespeedcf;
                //}


                LocalHeap lh(1000 * 1000);
                IntegrationRule ir (D==2?ET_TET:ET_TRIG, 0);
                MappedIntegrationPoint<3,3> mip(ir[0], ma->GetTrafo (ElementId(0), lh));

                this->gamma.SetSize(0);
                for(int i=0;i<ma->GetNE();i++) this->gamma.Append(Matrix<>(this->order-1));
                for(int ny=0;ny<=(this->order-2)*(D==2);ny++)
                {
                    for(int nx=0;nx<this->order-1;nx++)
                    {
                        double fac = (factorial(nx)*factorial(ny));
                        for(int ne=0;ne<ma->GetNE();ne++)
                        {
                            mip.Point() = ElCenter<2>(ElementId(ne));
                            this->gamma[ne](nx,ny) = localwavespeedcfx->Evaluate(mip)/fac;
                        }
                        localwavespeedcfx = localwavespeedcfx->Diff(MakeCoordinateCoefficientFunction(0).get(), make_shared<ConstantCoefficientFunction>(1) );
                    }
                    localwavespeedcf = localwavespeedcf->Diff(MakeCoordinateCoefficientFunction(1).get(), make_shared<ConstantCoefficientFunction>(1) );
                    localwavespeedcfx = localwavespeedcf;
                }
                cout << "finish" << endl;
            }
        }

        string GetClassName () const override { return "trefftz"; }

        void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

        FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

        size_t GetNDof () const override { return ndof; }

        static DocInfo GetDocu ();

        protected:

        template<int D>
        double Adiam(ElementId ei, double c) const;

        template<int D>
        double Adiam(ElementId ei, shared_ptr<CoefficientFunction> c) const;

        template<int D>
        Vec<D+1> ElCenter(ElementId ei) const;
    };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace(py::module m);
#endif // NGS_PYTHON

#endif
