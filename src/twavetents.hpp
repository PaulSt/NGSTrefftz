#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <tents.hpp>
#include "scalarmappedfe.hpp"
#include <unordered_map>

namespace ngcomp
{
  class TrefftzTents
  {
  public:
    TrefftzTents () { ; }
    virtual int dimensio () { return 0; }
    virtual void Propagate () { throw Exception ("TrefftzTents virtual!"); }
    virtual void SetInitial (shared_ptr<CoefficientFunction> init)
    {
      throw Exception ("TrefftzTents virtual!");
    }
    virtual void SetBoundaryCF (shared_ptr<CoefficientFunction> abddatum)
    {
      throw Exception ("TrefftzTents virtual!");
    }
  };

  template <int D> class TWaveTents : public TrefftzTents
  {
  protected:
    int order;
    shared_ptr<TentPitchedSlab> tps;
    shared_ptr<MeshAccess> ma;
    Vector<> wavespeed;
    shared_ptr<CoefficientFunction> wavespeedcf;
    Matrix<> wavefront;
    shared_ptr<CoefficientFunction> bddatum;
    int fosystem = 0;
    double timeshift = 0;
    int nbasis;
    const int nsimd = SIMD<double>::Size ();
    static constexpr ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);

    template <typename TFUNC>
    void
    CalcTentEl (int elnr, const Tent *tent, ScalarMappedElement<D + 1> &tel,
                TFUNC LocalWavespeed, SIMD_IntegrationRule &sir,
                LocalHeap &slh, SliceMatrix<> elmat, FlatVector<> elvec,
                SliceMatrix<SIMD<double>> simddshapes);

    void
    CalcTentBndEl (int surfel, const Tent *tent,
                   ScalarMappedElement<D + 1> &tel, SIMD_IntegrationRule &sir,
                   LocalHeap &slh, SliceMatrix<> elmat, FlatVector<> elvec);

    void CalcTentMacroEl (int fnr, const Array<int> &elnums,
                          std::unordered_map<int, int> &macroel,
                          const Tent *tent, ScalarMappedElement<D + 1> &tel,
                          SIMD_IntegrationRule &sir, LocalHeap &slh,
                          SliceMatrix<> elmat, FlatVector<> elvec);

    void
    CalcTentElEval (int elnr, const Tent *tent,
                    ScalarMappedElement<D + 1> &tel, SIMD_IntegrationRule &sir,
                    LocalHeap &slh, SliceVector<> sol,
                    SliceMatrix<SIMD<double>> simddshapes);

    Mat<D + 1, D + 1> TentFaceVerts (const Tent *tent, int elnr, int top);

    double TentFaceArea (Mat<D + 1, D + 1> v);

    Vec<D + 1> TentFaceNormal (Mat<D + 1, D + 1> v, int dir);

    template <typename T = double> void SwapIfGreater (T &a, T &b);

    double TentAdiam (const Tent *tent);

    inline void Solve (FlatMatrix<double> a, FlatVector<double> b);

    inline int MakeMacroEl (const Array<int> &tentel,
                            std::unordered_map<int, int> &macroel);

    void GetFacetSurfaceElement (shared_ptr<MeshAccess> ma, int fnr,
                                 Array<int> &selnums);

  public:
    TWaveTents (int aorder, shared_ptr<TentPitchedSlab> atps,
                double awavespeed)
        : order (aorder), tps (atps)
    {
      ma = atps->ma;
      nbasis
          = BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1);
      wavespeed.SetSize (1);
      wavespeed[0] = awavespeed;
      this->wavespeedcf
          = make_shared<ConstantCoefficientFunction> (awavespeed);
    }

    TWaveTents (int aorder, shared_ptr<TentPitchedSlab> atps,
                shared_ptr<CoefficientFunction> awavespeedcf)
        : order (aorder), tps (atps), wavespeedcf (awavespeedcf)
    {
      ma = atps->ma;
      nbasis
          = BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1);
      wavespeed.SetSize (ma->GetNE ());
      LocalHeap lh (1000 * 1000 * 1000);
      for (Ngs_Element el : ma->Elements (VOL))
        {
          ElementId ei = ElementId (el);
          // ELEMENT_TYPE eltype = ma->GetElType(ei);
          IntegrationRule ir (eltyp, 0);
          ElementTransformation &trafo = ma->GetTrafo (ei, lh);
          MappedIntegrationPoint<D, D> mip (ir[0], trafo);
          wavespeed[el.Nr ()] = awavespeedcf->Evaluate (mip);
        }
    }

    void Propagate () override;

    Matrix<>
    MakeWavefront (shared_ptr<CoefficientFunction> cf, double time = 0);

    Matrix<> GetWavefront () { return wavefront; }

    void SetInitial (shared_ptr<CoefficientFunction> init) override
    {
      wavefront = MakeWavefront (init);
      if (init->Dimension () == D + 1)
        {
          fosystem = 1;
          nbasis = BinCoeff (D + order, order)
                   + BinCoeff (D + order - 1, order - 1) - 1;
        }
    }

    void SetBoundaryCF (shared_ptr<CoefficientFunction> abddatum) override
    {
      bddatum = abddatum;
    }

    double Error (Matrix<> wavefront, Matrix<> wavefront_corr);

    double L2Error (Matrix<> wavefront, Matrix<> wavefront_corr);

    double Energy (Matrix<> wavefront);

    double MaxAdiam ();

    int LocalDofs () { return nbasis; }

    int GetOrder () { return order; }
    int GetSpaceDim () { return D; }
    shared_ptr<MeshAccess> GetInitmesh () { return ma; }
  };

  template <int D> class QTWaveTents : public TWaveTents<D>
  {
  private:
    Matrix<shared_ptr<CoefficientFunction>> GGder;
    Matrix<shared_ptr<CoefficientFunction>> BBder;
    double TentXdiam (const Tent *tent);
    const int nsimd = SIMD<double>::Size ();
    static constexpr ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);

    using TWaveTents<D>::Solve;
    using TWaveTents<D>::TentFaceVerts;

  public:
    QTWaveTents (int aorder, shared_ptr<TentPitchedSlab> atps,
                 shared_ptr<CoefficientFunction> awavespeedcf,
                 shared_ptr<CoefficientFunction> aBBcf)
        : TWaveTents<D> (aorder, atps, awavespeedcf)
    {
      this->nbasis = BinCoeff (D + this->order, this->order)
                     + BinCoeff (D + this->order - 1, this->order - 1);
      shared_ptr<CoefficientFunction> GGcf
          = make_shared<ConstantCoefficientFunction> (1)
            / (awavespeedcf * awavespeedcf);
      shared_ptr<CoefficientFunction> GGcfx
          = make_shared<ConstantCoefficientFunction> (1)
            / (awavespeedcf * awavespeedcf);
      GGder.SetSize (this->order - 1, (this->order - 2) * (D == 2) + 1);
      for (int ny = 0; ny <= (this->order - 2) * (D == 2); ny++)
        {
          for (int nx = 0; nx <= this->order - 2; nx++)
            {
              GGder (nx, ny) = GGcfx;
              GGcfx
                  = GGcfx->Diff (MakeCoordinateCoefficientFunction (0).get (),
                                 make_shared<ConstantCoefficientFunction> (1));
            }
          GGcf = GGcf->Diff (MakeCoordinateCoefficientFunction (1).get (),
                             make_shared<ConstantCoefficientFunction> (1));
          GGcfx = GGcf;
        }

      if (!aBBcf)
        aBBcf = make_shared<ConstantCoefficientFunction> (1);
      shared_ptr<CoefficientFunction> BBcf = aBBcf;
      shared_ptr<CoefficientFunction> BBcfx = aBBcf;
      BBder.SetSize (this->order, (this->order - 1) * (D == 2) + 1);
      for (int ny = 0; ny <= (this->order - 1) * (D == 2); ny++)
        {
          for (int nx = 0; nx <= this->order - 1; nx++)
            {
              BBder (nx, ny) = BBcfx;
              BBcfx
                  = BBcfx->Diff (MakeCoordinateCoefficientFunction (0).get (),
                                 make_shared<ConstantCoefficientFunction> (1));
            }
          BBcf = BBcf->Diff (MakeCoordinateCoefficientFunction (1).get (),
                             make_shared<ConstantCoefficientFunction> (1));
          BBcfx = BBcf;
        }
    }

    void Propagate ();
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTWaveTents (py::module m);
#endif // NGS_PYTHON

#endif
