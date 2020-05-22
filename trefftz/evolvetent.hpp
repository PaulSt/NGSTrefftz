#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
//#include <comp.hpp>    // provides FESpace, ...
#include <solve.hpp>
#include <h1lofe.hpp>
#include <fem.hpp>
#include "tents/tents.hpp"
#include "trefftzwavefe.hpp"
#include "trefftzgppwfe.hpp"


namespace ngcomp
{
    class TrefftzTents
    {
        private:

        public:
            TrefftzTents(){;}
            //virtual ~TrefftzTents() = default;
            virtual int dimensio(){return 0;}
    };

    template<int D>
    class WaveTents : public TrefftzTents
    {
        protected:
            int order;
            shared_ptr<MeshAccess> ma;
            Vector<> wavespeed;
            shared_ptr<CoefficientFunction> wavespeedcf;
            Matrix<> wavefront;
            shared_ptr<CoefficientFunction> bddatum;
            double timeshift = 0;

            void CalcTentEl(int elnr, Tent* tent, ScalarMappedElement<D+1> &tel,
                    SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes);

            void CalcTentBndEl(int surfel, Tent* tent, ScalarMappedElement<D+1> &tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec);

            void CalcTentMacroEl(int fnr, const Array<int> &elnums, std::unordered_map<int,int> &macroel, Tent* tent, TrefftzWaveFE<D> &tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec);

            void CalcTentElEval(int elnr, Tent* tent, ScalarMappedElement<D+1> &tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceVector<> sol, SliceMatrix<SIMD<double>> simddshapes);

            Mat<D+1,D+1> TentFaceVerts(Tent* tent, int elnr, int top);

            void TentDmat(Mat<D+1> &Dmat, Mat<D+1> v, int top, double wavespeed);

            double TentFaceArea( Mat<D+1,D+1> v );

            Vec<D+1> TentFaceNormal( Mat<D+1,D+1> v, int dir );

            template<typename T=double>
            void SwapIfGreater(T& a, T& b);

            double TentAdiam(Tent* tent);

            inline void LapackSolve(SliceMatrix<double> a, SliceVector<double> b);

            inline int MakeMacroEl(const Array<int> &tentel, std::unordered_map<int,int> &macroel);

        public:
            WaveTents(){;}

            WaveTents( int aorder, shared_ptr<MeshAccess> ama, double awavespeed, shared_ptr<CoefficientFunction> abddatum)
                : order(aorder), ma(ama), bddatum(abddatum)
            {
                wavespeed.SetSize(1);
                wavespeed[0]=awavespeed;
            }

            WaveTents( int aorder, shared_ptr<MeshAccess> ama, shared_ptr<CoefficientFunction> awavespeedcf, shared_ptr<CoefficientFunction> abddatum)
                : order(aorder), ma(ama), bddatum(abddatum), wavespeedcf(awavespeedcf)
            {
                wavespeed.SetSize(ama->GetNE());
                LocalHeap lh(1000 * 1000);
                for (Ngs_Element el : ama->Elements(VOL))
                {
                    ElementId ei = ElementId(el);
                    ELEMENT_TYPE eltype = ama->GetElType(ei);
                    IntegrationRule ir (eltype, 0);
                    ElementTransformation & trafo = ama->GetTrafo (ei, lh);
                    MappedIntegrationPoint<D,D> mip(ir[0], trafo);
                    wavespeed[el.Nr()] = awavespeedcf->Evaluate(mip);
                }
            }

            void EvolveTents(double dt);

            Matrix<> MakeWavefront( shared_ptr<CoefficientFunction> bddatum, double time);

            Matrix<> GetWavefront() {return wavefront;}
            void SetWavefront(shared_ptr<CoefficientFunction> bddatum, double time) { wavefront = MakeWavefront( bddatum, time); }
            //void SetWavefront(Matrix<> wf) { wavefront = wf; }

            double Error(Matrix<> wavefront, Matrix<> wavefront_corr);

            double L2Error(Matrix<> wavefront, Matrix<> wavefront_corr);

            double Energy(Matrix<> wavefront);

            double MaxAdiam(double dt);

            int LocalDofs(){
                TrefftzWaveFE<D> tel(order,wavespeed[0]);
                return tel.GetNDof();
            }

            int NrTents(double dt)
            {
                LocalHeap lh(1000 * 1000 * 100);
                TentPitchedSlab<D> tps = TentPitchedSlab<D>(ma);
                tps.PitchTents(dt, wavespeed[0]+1,lh);
                return tps.tents.Size();
            }

            int GetOrder(){return order;}
            int GetSpaceDim(){return D;}
            shared_ptr<MeshAccess> GetInitmesh(){return ma;}
    };


    template<int D>
    class GppwTents : public WaveTents<D>
    {
        private:
            //int order;
            //shared_ptr<MeshAccess> ma;
            //Vector<> wavespeed;
            //shared_ptr<CoefficientFunction> wavespeedcf;
            //Matrix<> wavefront;
            //shared_ptr<CoefficientFunction> bddatum;
            //double timeshift = 0;
            Array<Matrix<double>> gamma;

            using WaveTents<D>::TentAdiam;
            using WaveTents<D>::LapackSolve;
            using WaveTents<D>::TentFaceVerts;

        public:
            GppwTents( int aorder, shared_ptr<MeshAccess> ama, shared_ptr<CoefficientFunction> awavespeedcf, shared_ptr<CoefficientFunction> abddatum, shared_ptr<CoefficientFunction> x, shared_ptr<CoefficientFunction> y)
                : WaveTents<D>(aorder,ama,awavespeedcf,abddatum)
            {
                LocalHeap lh(1000 * 1000);
                const ELEMENT_TYPE eltyp = (D==3) ? ET_TET : ((D==2) ? ET_TRIG : ET_SEGM);
                const int nsimd = SIMD<double>::Size();
                SIMD_IntegrationRule sir(eltyp, this->order*2);

                IntegrationRule ir (eltyp, 0);
                for(int nv=0;nv<ama->GetNV();nv++)
                {
                    shared_ptr<CoefficientFunction> localwavespeedcf = awavespeedcf;
                    shared_ptr<CoefficientFunction> localwavespeedcfx = awavespeedcf;
                    MappedIntegrationPoint<D,D> mip(ir[0], ama->GetTrafo (ElementId(0), lh));
                    mip.Point() = ama->GetPoint<D>(nv);
                    Matrix<> b(this->order,this->order);
                    for(int nx=0;nx<this->order;nx++)
                    {
                        int ny = 0;
                        for(int ny=0;ny<this->order;ny++)
                        {
                            b(nx,ny) = localwavespeedcfx->Evaluate(mip);
                            localwavespeedcfx = localwavespeedcfx->Diff(y.get(), make_shared<ConstantCoefficientFunction>(1) );
                        }
                        localwavespeedcf = localwavespeedcf->Diff(x.get(), make_shared<ConstantCoefficientFunction>(1) );
                        localwavespeedcfx = localwavespeedcf;
                    }
                    this->gamma.Append(b);
                }
            }

            void EvolveTents(double dt);

            void CalcTentEl(int elnr, Tent* tent, ScalarMappedElement<D+1> &tel,
                    SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes);

            constexpr int factorial(int n)
            {
                return n>1 ? n * factorial(n-1) : 1;
            }
    };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m);
#endif // NGS_PYTHON

#endif
