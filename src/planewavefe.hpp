#ifndef FILE_PLANEWAVEELEMENT_HPP
#define FILE_PLANEWAVEELEMENT_HPP

#include <fem.hpp>
#include "ngsttd.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{
  template <int D> class PlaneWaveElement : public ScalarMappedElement<D>
  {
  private:
    Vec<D> GetDirection (int i) const;
    bool iscomplex = true;
    double elsize;
    double c;
    int conj;

  public:
    PlaneWaveElement (int andof, int aord, ELEMENT_TYPE aeltype,
                      Vec<D> ashift = 0, double aelsize = 1, double ac = 1.0,
                      int aconj = 1)
        : ScalarMappedElement<D> (andof, aord, Matrix<> (), aeltype, ashift,
                                  1.0),
          elsize (aelsize), c (ac), conj (aconj)
    {
      ;
    }

    using ScalarMappedElement<D>::CalcMappedDShape;
    using ScalarMappedElement<D>::Evaluate;
    using ScalarMappedElement<D>::EvaluateGrad;
    using ScalarMappedElement<D>::AddGradTrans;

    bool ComplexShapes () const override { return true; }

    Vec<D> EvaluateGrad (const BaseMappedIntegrationPoint &,
                         BareSliceVector<>) const override
    {
      throw Exception ("EvaluateGrad only complex for PW ");
    }
    void CalcShape (const BaseMappedIntegrationPoint &,
                    BareSliceVector<>) const override
    {
      throw Exception ("CalcShape only complex for PW ");
    }
    void CalcShape (const BaseMappedIntegrationRule &,
                    BareSliceMatrix<>) const override
    {
      throw Exception ("CalcShape only complex for PW ");
    }
    void CalcDShape (const BaseMappedIntegrationPoint &,
                     BareSliceMatrix<>) const override
    {
      throw Exception ("CalcDShape only complex for PW ");
    }
    void CalcDShape (const BaseMappedIntegrationRule &,
                     BareSliceMatrix<>) const override
    {
      throw Exception ("CalcDShape only complex for PW ");
    }
    void CalcMappedDDShape (const BaseMappedIntegrationPoint &,
                            BareSliceMatrix<>) const override
    {
      throw Exception ("Not implemented for PW ");
    }
    void CalcShape (const SIMD_BaseMappedIntegrationRule &,
                    BareSliceMatrix<SIMD<double>>) const override
    {
      throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
    }
    void CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                     BareSliceMatrix<SIMD<double>>) const override
    {
      throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
    }

    NGST_DLL virtual void
    Evaluate (const BaseMappedIntegrationRule &mir,
              BareSliceVector<Complex> coefs, FlatVector<Complex> vals) const;
    NGST_DLL virtual Complex
    EvaluateComplex (const BaseMappedIntegrationPoint &mip,
                     BareSliceVector<Complex> x) const;
    void CalcShape (const BaseMappedIntegrationPoint &mip,
                    BareSliceVector<Complex> shape) const override;
    void CalcDShape (const BaseMappedIntegrationPoint &mip,
                     BareSliceMatrix<Complex> dshape) const;

    void CalcDShape (const BaseMappedIntegrationRule &mir,
                     BareSliceMatrix<Complex> dshapes) const
    {
      for (size_t i = 0; i < mir.Size (); i++)
        CalcDShape (mir[i], dshapes.Cols (i * D, (i + 1) * D));
    }
    Vec<D> EvaluateGrad (const BaseMappedIntegrationPoint &,
                         BareSliceVector<Complex>) const
    {
      throw Exception ("OII");
    }

    INLINE const FlatMatrixFixWidth<D, Complex>
    GetDShapeComplex (const BaseMappedIntegrationPoint &mip,
                      LocalHeap &lh) const
    {
      FlatMatrixFixWidth<D, Complex> dshape (this->ndof, lh);
      CalcDShape (mip, dshape);
      return dshape;
    }
    virtual Vec<D, Complex>
    EvaluateGradComplex (const BaseMappedIntegrationPoint &ip,
                         BareSliceVector<Complex> x) const;
  };

  /// Identity
  template <int D, typename FEL = PlaneWaveElement<D>>
  class DiffOpMappedComplex : public DiffOp<DiffOpMappedComplex<D, FEL>>
  {
  public:
    enum
    {
      DIM = 1
    };
    enum
    {
      DIM_SPACE = D
    };
    enum
    {
      DIM_ELEMENT = D
    };
    enum
    {
      DIM_DMAT = 1
    };
    enum
    {
      DIFFORDER = 0
    };

    static string Name () { return "Id"; }
    static constexpr bool SUPPORT_PML = true;

    static const FEL &Cast (const FiniteElement &fel)
    {
      return static_cast<const FEL &> (fel);
    }

    static void
    GenerateMatrix (const FiniteElement &, const BaseMappedIntegrationPoint &,
                    BareSliceMatrix<double, ColMajor>, LocalHeap &)
    {
      throw Exception ("Not implemented for complex DiffOp");
    }
    static void
    GenerateMatrix (const FiniteElement &fel,
                    const BaseMappedIntegrationPoint &mip,
                    BareSliceMatrix<Complex, ColMajor> mat, LocalHeap &)
    {
      Cast (fel).CalcShape (mip, mat.Row (0));
    }

    template <typename MAT>
    static void GenerateMatrixIR (const FiniteElement &fel,
                                  const BaseMappedIntegrationRule &mir,
                                  MAT &mat, LocalHeap &)
    {
      Cast (fel).CalcShape (mir, Trans (mat));
    }

    template <typename MIP, class TVX, class TVY>
    static void
    Apply (const FiniteElement &, const MIP &, const TVX &, TVY &, LocalHeap &)
    {
      throw Exception ("Not implemented for complex DiffOp");
    }

    static void
    Apply (const FiniteElement &fel, const MappedIntegrationPoint<D, D> &mip,
           const BareSliceVector<Complex> &x, FlatVector<Complex> &y,
           LocalHeap &)
    {
      y (0) = Cast (fel).EvaluateComplex (mip, x);
    }

    template <class MIR, class TMY>
    static void ApplyIR (const FiniteElement &, const MIR &,
                         BareSliceVector<double>, TMY, LocalHeap &)
    {
      throw Exception ("Not implemented for complex DiffOp");
    }

    template <class MIR>
    static void
    ApplyIR (const FiniteElement &fel, const MIR &mir,
             BareSliceVector<Complex> x, SliceMatrix<Complex> y, LocalHeap &)
    {
      Cast (fel).Evaluate (mir, x,
                           FlatVector<Complex> (mir.Size (), &y (0, 0)));
    }

    template <typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FiniteElement &, const MIP &, const TVX &,
                            TVY &, LocalHeap &)
    {
      throw Exception ("Not implemented for complex DiffOp");
    }

    template <typename MIR, class TVX, class TVY>
    static void ApplyTransIR (const FiniteElement &fel, const MIR &mir,
                              const TVX &x, TVY &y, LocalHeap &lh)
    {
      y.Range (0, DIM * fel.GetNDof ()) = 0.0;
      for (size_t i = 0; i < mir.Size (); i++)
        {
          HeapReset hr (lh);
          ApplyTransAdd (fel, mir[i], x.Row (i), y, lh);
        }
    }

    template <typename MIP, class TVX, class TVY>
    static void ApplyTransAdd (const FiniteElement &fel, const MIP &mip,
                               const TVX &x, TVY &y, LocalHeap &lh)
    {
      HeapReset hr (lh);
      typedef typename TVX::TSCAL TSCAL;
      FlatMatrixFixHeight<DIM_DMAT, TSCAL> mat (DIM * fel.GetNDof (), lh);
      // Cast(fel).CalcShape (mip, mat.Row(0));
      GenerateMatrix (fel, mip, mat, lh);
      // y.Range (DIM * fel.GetNDof ()) += Trans (mat) * x;

      FlatVector<TSCAL> vec (DIM * fel.GetNDof (), lh);
      vec = Trans (mat) * x;
      y.Range (DIM * fel.GetNDof ()) += vec;
    }
  };

  /// Gradient operator of dimension D
  template <int D, typename FEL = PlaneWaveElement<D>>
  class DiffOpMappedGradientComplex
      : public DiffOp<DiffOpMappedGradientComplex<D, FEL>>
  {
  public:
    enum
    {
      DIM = 1
    };
    enum
    {
      DIM_SPACE = D
    };
    enum
    {
      DIM_ELEMENT = D
    };
    enum
    {
      DIM_DMAT = D
    };
    enum
    {
      DIFFORDER = 1
    };

    static string Name () { return "grad"; }
    static constexpr bool SUPPORT_PML = true;

    static const FEL &Cast (const FiniteElement &fel)
    {
      return static_cast<const FEL &> (fel);
    }

    template <typename SCALMIP, typename MAT>
    static void
    GenerateMatrix (const FiniteElement &fel,
                    const MappedIntegrationPoint<D, D, SCALMIP> &mip,
                    MAT &&mat, LocalHeap &)
    {
      Cast (fel).CalcDShape (mip, Trans (mat));
    }

    template <typename MIP, class TVY>
    static void
    Apply (const FiniteElement &fel, const MIP &mip,
           const BareSliceVector<Complex> &x, TVY &&y, LocalHeap &lh)
    {
      Vec<D, Complex> hv = Cast (fel).EvaluateGradComplex (mip, x);
      y = hv;
      HeapReset hr (lh);
    }

    template <typename MIP, class TVX, class TVY>
    static void Apply (const FiniteElement &, const MIP &, const TVX &, TVY &&,
                       LocalHeap &)
    {
      throw Exception ("Not implemented for complex DiffOp ");
    }

    // using DiffOp<DiffOpMappedGradient<D, FEL> >::ApplyTans;
    template <typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FiniteElement &fel, const MIP &mip,
                            const TVX &x, TVY &y, LocalHeap &lh)
    {
      typedef typename TVX::TSCAL TSCAL;
      Vec<D, TSCAL> vx = x;
      // auto hv = mip.GetJacobianInverse() * vx;
      y.Range (0, fel.GetNDof ())
          = Cast (fel).GetDShape (mip, lh) * vx; //* hv;
    }

    template <typename MIP, class TVX, class TVY>
    static void ApplyTransAdd (const FiniteElement &, const MIP &, const TVX &,
                               TVY &, LocalHeap &)
    {
      throw Exception ("Not implemented for complex DiffOp ");
    }
    // for complex shapes in lfi, symbolicintegrator.cpp T_CalcFacetVector ->
    // diffop_impl.hpp ApplyTrans,ApplyTransIR ->
    template <class MIP>
    static void ApplyTransAdd (const FiniteElement &fel, const MIP &mip,
                               FlatVector<Complex> x,
                               BareSliceVector<Complex> y, LocalHeap &lh)
    {
      FlatVector<Complex> vec (DIM * fel.GetNDof (), lh);
      vec = (Cast (fel).GetDShapeComplex (mip, lh)) * x;
      y.Range (DIM * fel.GetNDof ()) += vec;
    }

    template <typename MIR, class TVX, class TVY>
    static void ApplyTransIR (const FiniteElement &fel, const MIR &mir,
                              const TVX &x, TVY &y, LocalHeap &lh)
    {
      y.Range (0, DIM * fel.GetNDof ()) = 0.0;
      for (size_t i = 0; i < mir.Size (); i++)
        {
          HeapReset hr (lh);
          ApplyTransAdd (fel, mir[i], x.Row (i), y, lh);
        }
    }
  };

  template <typename DIFFOP>
  class T_DifferentialOperatorC : public T_DifferentialOperator<DIFFOP>
  {
    // copy past form diffop_impl to fix fast compile missing complex functions

    void
    Apply (const FiniteElement &bfel, const BaseMappedIntegrationRule &bmir,
           BareSliceVector<Complex> x, BareSliceMatrix<Complex> flux,
           LocalHeap &lh) const
    {
      auto fluxsize = flux.AddSize (bmir.Size (), DIFFOP::DIM_DMAT);
      const MappedIntegrationRule<DIFFOP::DIM_ELEMENT, DIFFOP::DIM_SPACE> &mir
          = static_cast<const MappedIntegrationRule<DIFFOP::DIM_ELEMENT,
                                                    DIFFOP::DIM_SPACE> &> (
              bmir);
      DIFFOP::ApplyIR (bfel, mir, x, fluxsize, lh);
    }

    void Apply (const FiniteElement &bfel,
                const SIMD_BaseMappedIntegrationRule &bmir,
                BareSliceVector<Complex> x,
                BareSliceMatrix<SIMD<Complex>> flux) const
    {
      DIFFOP::ApplySIMDIR (bfel, bmir, x, flux);
    }

    void ApplyTrans (const FiniteElement &bfel,
                     const BaseMappedIntegrationRule &bmir,
                     FlatMatrix<Complex> flux, BareSliceVector<Complex> x,
                     LocalHeap &lh) const
    {
      const MappedIntegrationRule<DIFFOP::DIM_ELEMENT, DIFFOP::DIM_SPACE> &mir
          = static_cast<const MappedIntegrationRule<DIFFOP::DIM_ELEMENT,
                                                    DIFFOP::DIM_SPACE> &> (
              bmir);
      DIFFOP::ApplyTransIR (bfel, mir, flux, x, lh);
    }
  };

}
#endif
