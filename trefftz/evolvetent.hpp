#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
//#include <comp.hpp>    // provides FESpace, ...
#include <solve.hpp>
#include <h1lofe.hpp>
#include "tents/tents.hpp"
#include "trefftzwavefe.hpp"
#include "trefftzgppwfe.hpp"

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
            Vector<> wavespeed;
            shared_ptr<CoefficientFunction> wavespeedcf;
            Matrix<> wavefront;
            shared_ptr<CoefficientFunction> bddatum;
            double timeshift = 0;

            void CalcTentEl(int elnr, Tent* tent, ScalarMappedElement<D+1> &tel,
                 SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes);

            void CalcTentBndEl(int surfel, Tent* tent, ScalarMappedElement<D+1> &tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec);

            void CalcTentElEval(int elnr, Tent* tent, ScalarMappedElement<D+1> &tel, SIMD_IntegrationRule &sir, LocalHeap &slh, SliceVector<> sol, SliceMatrix<SIMD<double>> simddshapes);

            Mat<D+1,D+1> TentFaceVerts(Tent* tent, int elnr, int top);

            void TentDmat(Mat<D+1> &Dmat, Mat<D+1> v, int top, double wavespeed);

            double TentFaceArea( Mat<D+1,D+1> v );

            Vec<D+1> TentFaceNormal( Mat<D+1,D+1> v, int dir );

            template<typename T=double>
            void SwapIfGreater(T& a, T& b);

            double TentAdiam(Tent* tent);

            inline void LapackSolve(SliceMatrix<double> a, SliceVector<double> b);

        public:
            WaveTents(){;}

            WaveTents( int aorder, shared_ptr<MeshAccess> ama, double awavespeed, shared_ptr<CoefficientFunction> abddatum)
                : order(aorder), ma(ama), bddatum(abddatum)
            {
                wavespeed.SetSize(1);
                wavespeed[0]=awavespeed;
            }

            //WaveTents( int aorder, shared_ptr<MeshAccess> ama, Vector<> awavespeed, shared_ptr<CoefficientFunction> abddatum)
                //: order(aorder), ma(ama), bddatum(abddatum), wavespeed(awavespeed)
            //{ ; }

            WaveTents( int aorder, shared_ptr<MeshAccess> ama, shared_ptr<CoefficientFunction> awavespeedcf, shared_ptr<CoefficientFunction> abddatum)
                : order(aorder), ma(ama), bddatum(abddatum), wavespeedcf(awavespeedcf)
            {
                wavespeed.SetSize(ma->GetNE());
                LocalHeap lh(1000 * 1000);
                for (Ngs_Element el : ma->Elements(VOL))
                {
                    ElementId ei = ElementId(el);
                    ELEMENT_TYPE eltype = ma->GetElType(ei);
                    IntegrationRule ir (eltype, 0);
                    ElementTransformation & trafo = ma->GetTrafo (ei, lh);
                    MappedIntegrationPoint<D,D> mip(ir[0], trafo);
                    wavespeed[el.Nr()] = wavespeedcf->Evaluate(mip);
                }
            }

            void EvolveTents(double dt);

            Matrix<> MakeWavefront( shared_ptr<CoefficientFunction> bddatum, double time);

            Matrix<> GetWavefront() {return wavefront;}
            void SetWavefront(shared_ptr<CoefficientFunction> bddatum, double time) { wavefront = MakeWavefront( bddatum, time); }
            //void SetWavefront(Matrix<> wf) { wavefront = wf; }

            double Error(Matrix<> wavefront, Matrix<> wavefront_corr);

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

            inline int MakeMacroEl(const Array<int> &tentel, std::unordered_map<int,int> &macroel)
            {
                // TODO fix if macro elements do not share faces
                int nrmacroel = 0;
                for(int i=0;i<tentel.Size();i++)
                {
                    int j=0;
                    while(wavespeed[tentel[i]]!=wavespeed[tentel[j]]) j++;
                    if(j==i)
                        macroel[tentel[i]] = nrmacroel++;
                    else
                        macroel[tentel[i]] = macroel[tentel[j]];
                }
                return nrmacroel;
            }
    };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m);
#endif // NGS_PYTHON

#endif
