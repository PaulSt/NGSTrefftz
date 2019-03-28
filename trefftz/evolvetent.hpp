#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
//#include <comp.hpp>    // provides FESpace, ...
#include <solve.hpp>
#include <h1lofe.hpp>
#include "tents/tents.hpp"
#include "trefftzwavefe.hpp"

namespace ngcomp
{

    class TrefftzTents {
        private:

        public:
            TrefftzTents(){;}
            //virtual ~TrefftzTents() = default;
            virtual int dimensio(){return 0;}

    };

    template<int D>
    class WaveTents : public TrefftzTents
    {
        private:
            int order;
            shared_ptr<MeshAccess> ma;
            double wavespeed;
            //Matrix<> wavefront;
            shared_ptr<CoefficientFunction> bddatum;
            double timeshift = 0;

            void CalcTentEl(int elnr, Tent* tent, TrefftzWaveFE<D+1> tel,
                    //SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec);
                 SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes, SliceMatrix<double> wavefront);

            void CalcTentBndEl(int surfel, Tent* tent, TrefftzWaveFE<D+1> tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec);

            //void CalcTentElEval(int elnr, Tent* tent, TrefftzWaveFE<D+1> tel, shared_ptr<MeshAccess> ma , SliceMatrix<> &wavefront, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceVector<> sol);
            void CalcTentElEval(int elnr, Tent* tent, TrefftzWaveFE<D+1> tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceVector<> sol, SliceMatrix<SIMD<double>> simddshapes, SliceMatrix<double> wavefront);

            Mat<D+1,D+1> TentFaceVerts(Tent* tent, int elnr, int top);

            void TentDmat(Mat<D+1> &Dmat, Mat<D+1> v, int top);

            double TentFaceArea( Mat<D+1,D+1> v );

            Vec<D+1> TentFaceNormal( Mat<D+1,D+1> v, int dir );

            template<typename T=double>
            void SwapIfGreater(T& a, T& b);

            double TentAdiam(Tent* tent);

            inline void LapackSolve(SliceMatrix<double> a, SliceVector<double> b);

        public:
            WaveTents(){;}

            WaveTents( int aorder, shared_ptr<MeshAccess> ama, double awavespeed, shared_ptr<CoefficientFunction> abddatum)
                : order(aorder), ma(ama), wavespeed(awavespeed), bddatum(abddatum)
            {
                //const ELEMENT_TYPE eltyp = (D==3) ? ET_TET : ((D==2) ? ET_TRIG : ET_SEGM );
                //IntegrationRule ir(eltyp, order*2);
                //int nsimd = SIMD<double>::Size();
                //int snip = ir.Size() + (ir.Size()%nsimd==0?0:nsimd-ir.Size()%nsimd);
                //wavefront.SetSize(ama->GetNE(),snip * (D+2));
                //wavefront = MakeWavefront(abddatum, 0);
                //cout << wavefront;
            }

            Matrix<> EvolveTents(double dt, Matrix<double> wavefront);

            Matrix<> MakeWavefront( shared_ptr<CoefficientFunction> bddatum, double time);

            double Error(Matrix<> wavefront, Matrix<> wavefront_corr);

            double Energy(Matrix<> wavefront);

            int dimensio(){return D;}
    };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m);
#endif // NGS_PYTHON

#endif
