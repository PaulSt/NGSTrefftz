#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
//#include <comp.hpp>    // provides FESpace, ...
#include <solve.hpp>
#include <h1lofe.hpp>
#include "tents/tents.hpp"
#include "trefftzwavefe.hpp"

namespace ngcomp
{

  class TrefftzTents
  {
  private:
  public:
    TrefftzTents () { ; }
    // virtual ~TrefftzTents() = default;
    virtual int dimensio () { return 0; }
  };

  template <int D> class WaveTents : public TrefftzTents
  {
  private:
    int order;
    shared_ptr<MeshAccess> ma;
    double wavespeed;
    // Matrix<> wavefront;
    shared_ptr<CoefficientFunction> bddatum;
    double timeshift = 0;

    void
    CalcTentEl (int elnr, Tent *tent, TrefftzWaveFE<D + 1> tel,
                SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat,
                SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes,
                Matrix<> &wavefront);

    void CalcTentBndEl (int surfel, Tent *tent, TrefftzWaveFE<D + 1> tel,
                        SIMD_IntegrationRule &sir, LocalHeap &slh,
                        SliceMatrix<> elmat, SliceVector<> elvec);

    void
    CalcTentElEval (int elnr, Tent *tent, TrefftzWaveFE<D + 1> tel,
                    SIMD_IntegrationRule &sir, LocalHeap &slh,
                    SliceVector<> sol, SliceMatrix<SIMD<double>> simddshapes,
                    Matrix<> &wavefront);

    Mat<D + 1, D + 1> TentFaceVerts (Tent *tent, int elnr, int top);

    void TentDmat (Mat<D + 1> &Dmat, Mat<D + 1> v, int top);

    double TentFaceArea (Mat<D + 1, D + 1> v);

    Vec<D + 1> TentFaceNormal (Mat<D + 1, D + 1> v, int dir);

    template <typename T = double> void SwapIfGreater (T &a, T &b);

    double TentAdiam (Tent *tent);

    inline void LapackSolve (SliceMatrix<double> a, SliceVector<double> b);

  public:
    WaveTents () { ; }

    WaveTents (int aorder, shared_ptr<MeshAccess> ama, double awavespeed,
               shared_ptr<CoefficientFunction> abddatum)
        : order (aorder), ma (ama), wavespeed (awavespeed), bddatum (abddatum)
    {
      ;
    }

    void EvolveTents (double dt, Matrix<> &wavefront);

    Matrix<>
    MakeWavefront (shared_ptr<CoefficientFunction> bddatum, double time);

    double Error (Matrix<> wavefront, Matrix<> wavefront_corr);

    double Energy (Matrix<> wavefront);

    double MaxAdiam (double dt);

    int LocalDofs ()
    {
      TrefftzWaveFE<D + 1> tel (order, wavespeed);
      return tel.GetNBasis ();
    }

    int NrTents (double dt)
    {
      TentPitchedSlab<D> tps = TentPitchedSlab<D> (ma);
      tps.PitchTents (dt, wavespeed + 1);
      return tps.tents.Size ();
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m);
#endif // NGS_PYTHON

#endif
