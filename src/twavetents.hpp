#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <tents.hpp>
#include "scalarmappedfe.hpp"
#include "trefftzfespace.hpp"
#include <unordered_map>

namespace ngfem
{
  template <int DIM_ELEMENT, int DIM_SPACE>
  class SIMD_STMappedIntegrationRule : public SIMD_BaseMappedIntegrationRule
  {
    FlatArray<SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>> mips;

  public:
    SIMD_STMappedIntegrationRule (const SIMD_IntegrationRule &ir,
                                  const ElementTransformation &aeltrans,
                                  Allocator &lh)
        : SIMD_BaseMappedIntegrationRule (ir, aeltrans), mips (ir.Size (), lh)
    {
      throw Exception ("Not implemented for sstmip");
    }

    SIMD_STMappedIntegrationRule (const SIMD_IntegrationRule &ir,
                                  const ElementTransformation &aeltrans, int,
                                  Allocator &lh)
        : SIMD_BaseMappedIntegrationRule (ir, aeltrans), mips (ir.Size (), lh)
    {
      dim_element = DIM_ELEMENT;
      dim_space = DIM_SPACE;
      baseip = (char *)(void *)(SIMD<BaseMappedIntegrationPoint> *)(&mips[0]);
      incr = sizeof (SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>);

      for (size_t i = 0; i < ir.Size (); i++)
        new (&mips[i]) SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>> (
            ir[i], eltrans, -1);

      new (&points) BareSliceMatrix<SIMD<double>> (
          mips.Size (), DIM_SPACE,
          sizeof (SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)
              / sizeof (SIMD<double>),
          &mips[0].Point () (0));

      new (&normals) BareSliceMatrix<SIMD<double>> (
          mips.Size (), DIM_SPACE,
          sizeof (SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)
              / sizeof (SIMD<double>),
          &mips[0].NV () (0));
    }

    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE, int) override
    {
      throw Exception ("Not implemented for sstmip");
    }
    SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>> &
    operator[] (size_t i) const
    {
      return mips[i];
    }
    virtual void Print (ostream &ost) const override;

    virtual void
    TransformGradient (BareSliceMatrix<SIMD<double>>) const override
    {
      throw Exception ("Not implemented for sstmip");
    }
    virtual void
    TransformGradientTrans (BareSliceMatrix<SIMD<double>>) const override
    {
      throw Exception ("Not implemented for sstmip");
    }
  };
}

namespace ngcomp
{

  class TrefftzTents
  {
  public:
    TrefftzTents () { ; }
    virtual ~TrefftzTents () = default;
    virtual int dimensio () { return 0; }
    virtual void Propagate () { throw Exception ("TrefftzTents virtual!"); }
    virtual void SetInitial (shared_ptr<CoefficientFunction>)
    {
      throw Exception ("TrefftzTents virtual!");
    }
    virtual void SetBoundaryCF (shared_ptr<CoefficientFunction>)
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
    const size_t nsimd = SIMD<double>::Size ();
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
    QTWaveBasis<D + 1> basis;
    // Matrix<shared_ptr<CoefficientFunction>> GGder;
    // Matrix<shared_ptr<CoefficientFunction>> BBder;
    double TentXdiam (const Tent *tent);
    const size_t nsimd = SIMD<double>::Size ();
    static constexpr ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);

    using TWaveTents<D>::Solve;
    using TWaveTents<D>::TentFaceVerts;

  public:
    QTWaveTents (int aorder, shared_ptr<TentPitchedSlab> atps,
                 shared_ptr<CoefficientFunction> awavespeedcf,
                 shared_ptr<CoefficientFunction> aBBcf)
        : TWaveTents<D> (aorder, atps, awavespeedcf),
          basis (aorder, awavespeedcf, aBBcf)
    {
      this->nbasis = BinCoeff (D + this->order, this->order)
                     + BinCoeff (D + this->order - 1, this->order - 1);
    }

    void Propagate ();
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTWaveTents (py::module m);
#endif // NGS_PYTHON

#endif
