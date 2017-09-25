#ifndef FILE_DIFFOPMAPPED
#define FILE_DIFFOPMAPPED

#include "TrefftzElement.hpp"

namespace ngfem
{

  /*
     realizations of bdb integrators for many equations.
     The differential operators provide the B-matrix,
     the DMatOps provide the coefficient tensors
  */

  /// Identity
  template <int D, int order, typename FEL = TrefftzElement<D, order>>
  class DiffOpMapped : public DiffOp<DiffOpMapped<D, order, FEL>>
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

    static string Name () { return "Trefftz"; }
    static constexpr bool SUPPORT_PML = true;

    static const FEL &Cast (const FiniteElement &fel)
    {
      return static_cast<const FEL &> (fel);
    }

    static void GenerateMatrix (const FiniteElement &fel,
                                const BaseMappedIntegrationPoint &mip,
                                FlatMatrixFixHeight<1> &mat, LocalHeap &lh)
    {
      Cast (fel).CalcShape (
          mip, mat.Row (0)); // FlatVector<> (fel.GetNDof(), &mat(0,0)));
    }

    static void
    GenerateMatrix (const FiniteElement &fel,
                    const BaseMappedIntegrationPoint &mip,
                    SliceMatrix<double, ColMajor> mat, LocalHeap &lh)
    {
      Cast (fel).CalcShape (mip, mat.Row (0));
    }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement &fel, const MIP &mip,
                                MAT &&mat, LocalHeap &lh)
    {
      HeapReset hr (lh);
      mat.Row (0) = Cast (fel).GetShape (mip, lh);
    }
    /*
                    // using DiffOp<DiffOpId<D, FEL> > :: GenerateMatrixIR;
                    template <typename MAT>
                    static void GenerateMatrixIR (const FiniteElement & fel,
                                                                                                                                            const BaseMappedIntegrationRule & mir,
                                                                                                                                            MAT & mat, LocalHeap & lh)
                    {
                            Cast(fel).CalcShape (mir.IR(), Trans(mat));
                    }

                    static void GenerateMatrixSIMDIR (const FiniteElement &
       fel, const SIMD_BaseMappedIntegrationRule & mir,
                                                                                                                                                            BareSliceMatrix<SIMD<double>> mat)
                    {
                            Cast(fel).CalcShape (mir.IR(), mat);
                    }

                    template <typename MIP, class TVX, class TVY>
                    static void Apply (const FiniteElement & fel, const MIP &
       mip, const TVX & x, TVY & y, LocalHeap & lh)
                    {
                            HeapReset hr(lh);
                            y = Trans (Cast(fel).GetShape (mip.IP(), lh)) * x;
                    }

                    static void Apply (const FiniteElement & fel, const
       MappedIntegrationPoint<D,D> & mip, const FlatVector<double> & x,
       FlatVector<double> & y, LocalHeap & lh)
                    {
                            y(0) = Cast(fel).Evaluate(mip.IP(), x);
                    }


                    // using DiffOp<DiffOpId<D, FEL> >::ApplyIR;

                    template <class MIR, class TMY>
                    static void ApplyIR (const FiniteElement & fel, const MIR &
       mir, FlatVector<double> x, TMY y, LocalHeap & lh)
                    {
                            Cast(fel).Evaluate (mir.IR(), x, FlatVector<>
       (mir.Size(), &y(0,0)));
                    }

                    template <class MIR>
                    static void ApplyIR (const FiniteElement & fel, const MIR &
       mir, FlatVector<Complex> x, FlatMatrix<Complex> y, LocalHeap & lh)
                    {
                            Cast(fel).Evaluate (mir.IR(),
                                                                                                            SliceMatrix<double> (fel.GetNDof(), 2, 2, reinterpret_cast<double*> (&x(0))),
                                                                                                            SliceMatrix<double> (mir.Size(), 2, 2, reinterpret_cast<double*> (&y(0))));
                    }

                    // using ApplySIMDIR;
                    using DiffOp<DiffOpId<D, FEL> >::ApplySIMDIR;
                    static void ApplySIMDIR (const FiniteElement & fel, const
       SIMD_BaseMappedIntegrationRule & mir, BareSliceVector<double> x,
       BareSliceMatrix<SIMD<double>> y)
                    {
                            Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
                    }



                    template <typename MIP, class TVX, class TVY>
                    static void ApplyTrans (const FiniteElement & fel, const
       MIP & mip, const TVX & x, TVY & y, LocalHeap & lh)
                    {
                            HeapReset hr(lh);
                            y = Cast(fel).GetShape (mip.IP(), lh) * x;
                    }


                    // using DiffOp<DiffOpId<D, FEL> >::ApplyTransIR;
                    template <class MIR>
                    static void ApplyTransIR (const FiniteElement & fel,
                                                    const MIR & mir,
                                                    FlatMatrix<double> x,
       FlatVector<double> y, LocalHeap & lh)
                    {
                            Cast(fel).EvaluateTrans (mir.IR(), FlatVector<>
       (mir.Size(), &x(0,0)), y);
                    }

                    template <class MIR>
                    static void ApplyTransIR (const FiniteElement & fel,
                                                    const MIR & mir,
                                                    FlatMatrix<Complex> x,
       FlatVector<Complex> y, LocalHeap & lh)
                    {
                            DiffOp<DiffOpId<D, FEL> > :: ApplyTransIR (fel,
       mir, x, y, lh);
                    }

                    using DiffOp<DiffOpId<D, FEL> >::AddTransSIMDIR;
                    static void AddTransSIMDIR (const FiniteElement & fel,
       const SIMD_BaseMappedIntegrationRule & mir,
                                                                                                                                    BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
                    {
                            Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
                    }
            */
  };

  /// Identity on boundary
  template <int D, int order, typename FEL = TrefftzElement<D, order>>
  class DiffOpMappedBoundary
      : public DiffOp<DiffOpMappedBoundary<D, order, FEL>>
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
      DIM_ELEMENT = D - 1
    };
    enum
    {
      DIM_DMAT = 1
    };
    enum
    {
      DIFFORDER = 0
    };

    static const FEL &Cast (const FiniteElement &fel)
    {
      return static_cast<const FEL &> (fel);
    }

    template <typename AFEL, typename MIP, typename MAT>
    static void
    GenerateMatrix (const AFEL &fel, const MIP &mip, MAT &mat, LocalHeap &lh)
    {
      HeapReset hr (lh);
      mat.Row (0) = static_cast<const FEL &> (fel).GetShape (mip, lh);
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL &fel, const MIP &mip, const TVX &x, TVY &y,
                       LocalHeap &lh)
    {
      HeapReset hr (lh);
      y = Trans (static_cast<const FEL &> (fel).GetShape (mip, lh)) * x;
      // y(0) = InnerProduct (x, static_cast<const FEL&>(fel).GetShape
      // (mip.IP(), lh));
    }

    static void
    Apply (const ScalarFiniteElement<D - 1> &fel,
           const MappedIntegrationPoint<D - 1, D> &mip,
           const FlatVector<double> &x, FlatVector<double> &y, LocalHeap &lh)
    {
      y (0) = static_cast<const FEL &> (fel).Evaluate (mip, x);
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const AFEL &fel, const MIP &mip, const TVX &x,
                            TVY &y, LocalHeap &lh)
    {
      HeapReset hr (lh);
      y = static_cast<const FEL &> (fel).GetShape (mip, lh) * x;
    }

    /*
        // using DiffOp<DiffOpIdBoundary<D, FEL> >::ApplyTransIR;
        template <class MIR>
        static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                                  FlatMatrix<double> x, FlatVector<double> y,
                                  LocalHeap & lh)
        {
          // static Timer t("applytransir - bnd");
          // RegionTimer reg(t);

          static_cast<const FEL&>(fel).
            EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
        }

        template <class MIR>
        static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                                  FlatMatrix<Complex> x, FlatVector<Complex> y,
                                  LocalHeap & lh)
        {
          DiffOp<DiffOpIdBoundary<D, FEL> > :: ApplyTransIR (fel, mir, x, y,
       lh);
          // static_cast<const FEL&>(fel).
          // EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
        }

        using DiffOp<DiffOpIdBoundary<D, FEL> >::ApplySIMDIR;
        static void ApplySIMDIR (const FiniteElement & fel, const
       SIMD_BaseMappedIntegrationRule & mir, BareSliceVector<double> x,
       BareSliceMatrix<SIMD<double>> y)
        {
          Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
        }

        using DiffOp<DiffOpIdBoundary<D, FEL> >::AddTransSIMDIR;
        static void AddTransSIMDIR (const FiniteElement & fel, const
       SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> y,
       BareSliceVector<double> x)
        {
          Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
        }
    */
  };

}
#endif
