#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
//#include <comp.hpp>    // provides FESpace, ...
#include <solve.hpp>
#include <h1lofe.hpp>
#include <fem.hpp>
#include "../tents/tents.hpp"
#include "trefftzwavefe.hpp"
#include "qtrefftzwavefe.hpp"

namespace ngcomp
{

  static int addtentslope = 3;

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
  protected:
    int order;
    shared_ptr<MeshAccess> ma;
    Vector<> wavespeed;
    shared_ptr<CoefficientFunction> wavespeedcf;
    Matrix<> wavefront;
    shared_ptr<CoefficientFunction> bddatum;
    double timeshift = 0;

    void
    CalcTentEl (int elnr, Tent *tent, ScalarMappedElement<D + 1> &tel,
                SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat,
                SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes);

    void
    CalcTentBndEl (int surfel, Tent *tent, ScalarMappedElement<D + 1> &tel,
                   SIMD_IntegrationRule &sir, LocalHeap &slh,
                   SliceMatrix<> elmat, SliceVector<> elvec);

    void
    CalcTentMacroEl (int fnr, const Array<int> &elnums,
                     std::unordered_map<int, int> &macroel, Tent *tent,
                     TrefftzWaveFE<D> &tel, SIMD_IntegrationRule &sir,
                     LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec);

    void
    CalcTentElEval (int elnr, Tent *tent, ScalarMappedElement<D + 1> &tel,
                    SIMD_IntegrationRule &sir, LocalHeap &slh,
                    SliceVector<> sol, SliceMatrix<SIMD<double>> simddshapes);

    Mat<D + 1, D + 1> TentFaceVerts (Tent *tent, int elnr, int top);

    void TentDmat (Mat<D + 1> &Dmat, Mat<D + 1> v, int top, double wavespeed);

    double TentFaceArea (Mat<D + 1, D + 1> v);

    Vec<D + 1> TentFaceNormal (Mat<D + 1, D + 1> v, int dir);

    template <typename T = double> void SwapIfGreater (T &a, T &b);

    double TentAdiam (Tent *tent);

    inline void LapackSolve (SliceMatrix<double> a, SliceVector<double> b);

    inline int MakeMacroEl (const Array<int> &tentel,
                            std::unordered_map<int, int> &macroel);

  public:
    WaveTents () { ; }

    WaveTents (int aorder, shared_ptr<MeshAccess> ama, double awavespeed,
               shared_ptr<CoefficientFunction> abddatum)
        : order (aorder), ma (ama), bddatum (abddatum)
    {
      wavespeed.SetSize (1);
      wavespeed[0] = awavespeed;
      this->wavespeedcf
          = make_shared<ConstantCoefficientFunction> (awavespeed);
    }

    WaveTents (int aorder, shared_ptr<MeshAccess> ama,
               shared_ptr<CoefficientFunction> awavespeedcf,
               shared_ptr<CoefficientFunction> abddatum)
        : order (aorder), ma (ama), bddatum (abddatum),
          wavespeedcf (awavespeedcf)
    {
      wavespeed.SetSize (ama->GetNE ());
      LocalHeap lh (1000 * 1000 * 1000);
      for (Ngs_Element el : ama->Elements (VOL))
        {
          ElementId ei = ElementId (el);
          ELEMENT_TYPE eltype = ama->GetElType (ei);
          IntegrationRule ir (eltype, 0);
          ElementTransformation &trafo = ama->GetTrafo (ei, lh);
          MappedIntegrationPoint<D, D> mip (ir[0], trafo);
          wavespeed[el.Nr ()] = awavespeedcf->Evaluate (mip);
        }
    }

    void EvolveTents (double dt);

    Matrix<>
    MakeWavefront (shared_ptr<CoefficientFunction> bddatum, double time);

    Matrix<> GetWavefront () { return wavefront; }
    void SetWavefront (shared_ptr<CoefficientFunction> bddatum, double time)
    {
      wavefront = MakeWavefront (bddatum, time);
    }
    // void SetWavefront(Matrix<> wf) { wavefront = wf; }

    double Error (Matrix<> wavefront, Matrix<> wavefront_corr);

    double L2Error (Matrix<> wavefront, Matrix<> wavefront_corr);

    double Energy (Matrix<> wavefront);

    double MaxAdiam (double dt);

    int LocalDofs ()
    {
      TrefftzWaveFE<D> tel (order, wavespeed[0]);
      return tel.GetNDof ();
    }

    int NrTents (double dt)
    {
      LocalHeap lh (1000 * 1000 * 1000);
      TentPitchedSlab<D> tps = TentPitchedSlab<D> (ma);
      tps.PitchTents (
          dt,
          this->wavespeedcf
              + make_shared<ConstantCoefficientFunction> (addtentslope),
          lh);
      return tps.tents.Size ();
    }

    int GetOrder () { return order; }
    int GetSpaceDim () { return D; }
    shared_ptr<MeshAccess> GetInitmesh () { return ma; }
  };

  template <int D> class QTWaveTents : public WaveTents<D>
  {
  private:
    // int order;
    // shared_ptr<MeshAccess> ma;
    // Vector<> wavespeed;
    // shared_ptr<CoefficientFunction> wavespeedcf;
    // Matrix<> wavefront;
    // shared_ptr<CoefficientFunction> bddatum;
    // double timeshift = 0;
    Array<Matrix<double>> gamma;

    using WaveTents<D>::TentAdiam;
    using WaveTents<D>::LapackSolve;
    using WaveTents<D>::TentFaceVerts;

  public:
    QTWaveTents (int aorder, shared_ptr<MeshAccess> ama,
                 shared_ptr<CoefficientFunction> awavespeedcf,
                 shared_ptr<CoefficientFunction> abddatum)
        : WaveTents<D> (aorder, ama, awavespeedcf, abddatum)
    {
      LocalHeap lh (1000 * 1000 * 1000);
      const ELEMENT_TYPE eltyp
          = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
      const int nsimd = SIMD<double>::Size ();
      SIMD_IntegrationRule sir (eltyp, this->order * 2);

      IntegrationRule ir (eltyp, 0);
      shared_ptr<CoefficientFunction> localwavespeedcf
          = make_shared<ConstantCoefficientFunction> (1)
            / (awavespeedcf * awavespeedcf);
      shared_ptr<CoefficientFunction> localwavespeedcfx
          = make_shared<ConstantCoefficientFunction> (1)
            / (awavespeedcf * awavespeedcf);

      // this->gamma.SetSize(ama->GetNV());
      // for(auto &m : this->gamma)
      //{
      // Matrix<> b(this->order,(this->order-1)*(D==2)+1);
      // m = b;
      // }

      this->gamma.SetSize (0);
      for (int i = 0; i < ama->GetNV (); i++)
        this->gamma.Append (Matrix<> (this->order - 1));
      MappedIntegrationPoint<D, D> mip (ir[0],
                                        ama->GetTrafo (ElementId (0), lh));
      for (int ny = 0; ny <= (this->order - 2) * (D == 2); ny++)
        {
          for (int nx = 0; nx < this->order - 1; nx++)
            {
              double fac = (factorial (nx) * factorial (ny));
              for (int nv = 0; nv < ama->GetNV (); nv++)
                {
                  mip.Point () = ama->GetPoint<D> (nv);
                  this->gamma[nv](nx, ny)
                      = localwavespeedcfx->Evaluate (mip) / fac;
                }
              localwavespeedcfx = localwavespeedcfx->Diff (
                  MakeCoordinateCoefficientFunction (0).get (),
                  make_shared<ConstantCoefficientFunction> (1));
            }
          localwavespeedcf = localwavespeedcf->Diff (
              MakeCoordinateCoefficientFunction (1).get (),
              make_shared<ConstantCoefficientFunction> (1));
          localwavespeedcfx = localwavespeedcf;
        }
    }

    QTWaveTents (int aorder, shared_ptr<MeshAccess> ama,
                 shared_ptr<CoefficientFunction> awavespeedcf,
                 shared_ptr<CoefficientFunction> abddatum,
                 vector<shared_ptr<CoefficientFunction>> taylorcf)
        : WaveTents<D> (aorder, ama, awavespeedcf, abddatum)
    {
      LocalHeap lh (1000 * 1000 * 1000);
      const ELEMENT_TYPE eltyp
          = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
      const int nsimd = SIMD<double>::Size ();
      SIMD_IntegrationRule sir (eltyp, this->order * 2);

      IntegrationRule ir (eltyp, 0);

      cout << "start" << ama->GetNV () << endl;
      for (int i = 0; i < ama->GetNV (); i++)
        this->gamma.Append (Matrix<> (this->order, this->order));
      MappedIntegrationPoint<D, D> mip (ir[0],
                                        ama->GetTrafo (ElementId (0), lh));
      for (int nx, count = 0; nx < this->order; nx++)
        {
          for (int ny = 0; ny <= (this->order - 1) * (D == 2); ny++)
            {
              double fac = (factorial (nx) * factorial (ny));
              for (int nv = 0; nv < ama->GetNV (); nv++)
                {
                  mip.Point () = ama->GetPoint<D> (nv);
                  this->gamma[nv](nx, ny)
                      = taylorcf[count]->Evaluate (mip) / fac;
                }
              count++;
            }
        }
      cout << "finish" << endl;
    }

    void EvolveTents (double dt);

    void
    CalcTentEl (int elnr, Tent *tent, ScalarMappedElement<D + 1> &tel,
                SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat,
                SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes);
  };

}
// computeone: ex [master !?*]$ cat ../tex/numerics/tentquadgppw1th1.csv
// hnr,h,p,error,cond,ndof,time,rate
// 1,0.5,3,0.00022520811379308172,1,63,0.0015499591827392578,12.11645357380905
// 2,0.25,3,4.0857230174787716e-05,1,189,0.0031397342681884766,2.462595497548919
// 3,0.125,3,5.026517792686095e-06,1,630,0.009650468826293945,3.02296020389384
// 4,0.0625,3,5.483089995654643e-07,1,2331,0.029244661331176758,3.1964982357300626
// 5,0.03125,3,8.433908219043716e-08,1,8876,0.07805562019348145,2.7007159263637415
// 6,0.015625,3,9.871428843411423e-09,1,34762,0.2570507526397705,3.0948704934474556
// 1,0.5,4,5.5720909108324206e-05,1,81,0.0010800361633300781,14.131421678174396
// 2,0.25,4,4.827484663381129e-06,1,243,0.003147125244140625,3.5288752153190917
// 3,0.125,4,3.1949555323869466e-07,1,810,0.010710716247558594,3.917403918935851
// 4,0.0625,4,1.4960830249816737e-08,1,2997,0.04068136215209961,4.416533707001708
// 5,0.03125,4,6.622642197540464e-10,1,11412,0.15599751472473145,4.49763951380215
// 6,0.015625,4,2.986814263906881e-11,1,44694,0.6007552146911621,4.470727484587121
// computeone: ex [master !?*]$ cat ../tex/numerics/tentquadgppw2th1.csv
// hnr,h,p,error,cond,ndof,time,rate
// 1,0.5,3,0.00347547742990016,1,576,0.013589620590209961,8.168573108714414
// 2,0.25,3,0.00037628669352707446,1,3136,0.09395980834960938,3.207306997963732
// 3,0.125,3,5.9560518307497784e-05,1,19824,0.6843569278717041,2.6594040591643235
// 4,0.0625,3,9.690693262687142e-06,1,148576,5.7606916427612305,2.6196845251618393

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m);
#endif // NGS_PYTHON

#endif
