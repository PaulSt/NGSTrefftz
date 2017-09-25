#ifndef FILE_MAPPEDELEMENT
#define FILE_MAPPEDELEMENT

#include <fem.hpp>

namespace ngfem
{

  class MappedElement : public FiniteElement
  {
  public:
    // using FiniteElement::FiniteElement;

    INLINE MappedElement () { ; }
    INLINE MappedElement (int andof, int aorder)
        : FiniteElement (andof, aorder)
    {
      ;
    }

    /// compute shape
    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationPoint &mip,
               BareSliceVector<> shape) const = 0;

    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationPoint &mip,
               BareSliceVector<Complex> shape) const;

    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER virtual void
    CalcDShape (const BaseMappedIntegrationPoint &mip,
                SliceMatrix<> dshape) const = 0;

    /**
       returns shape functions in point ip.
    */
    INLINE FlatVector<>
    GetShape (const BaseMappedIntegrationPoint &mip, LocalHeap &lh) const
    {
      FlatVector<> shape (ndof, lh);
      CalcShape (mip, shape);
      return shape;
    }

    /// compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationRule &mir,
               SliceMatrix<> shape) const;

    /// compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER virtual void
    CalcShape (const SIMD_BaseMappedIntegrationRule &mir,
               BareSliceMatrix<SIMD<double>> shape) const;

    // rows dim*ndof, cols .. nip
    // rows:  phi0/dx, phi0/dy, phi0/dz, phi1/dx ...
    /*
HD NGS_DLL_HEADER
virtual void CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir,
                       BareSliceMatrix<SIMD<double>> dshapes) const;

                                                                                                                             */
    /**
       Evaluates function in integration point ip.
       Vector x provides coefficient vector.
     */
    HD NGS_DLL_HEADER virtual double
    Evaluate (const BaseMappedIntegrationPoint &mip,
              BareSliceVector<> x) const;

    /*
Evaluate function in points of integrationrule ir.
Vector x provides coefficient vector.
*/
    // HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir,
    // BareSliceVector<> coefs, FlatVector<> values) const; HD NGS_DLL_HEADER
    // virtual void Evaluate (const SIMD_IntegrationRule & ir,
    // BareSliceVector<> coefs, BareVector<SIMD<double>> values) const; HD
    // NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir,
    // SliceMatrix<> coefs, BareSliceMatrix<SIMD<double>> values) const;

    /**
       Each column a vector ...
     */
    // HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir,
    // SliceMatrix<> coefs, SliceMatrix<> values) const;

    /**
       Evaluate function in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */
    /*
        HD NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule &
       ir, FlatVector<> values, BareSliceVector<> coefs) const; HD
       NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir,
       BareVector<SIMD<double>> values, BareSliceVector<> coefs) const; HD
       NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir,
       BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const;

        HD NGS_DLL_HEADER virtual void EvaluateGrad (const
       SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs,
       BareSliceMatrix<SIMD<double>> values) const;
        // needed for ALE-trafo
        HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_IntegrationRule
       & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values)
       const; HD NGS_DLL_HEADER virtual void AddGradTrans (const
       SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>>
       values, BareSliceVector<> coefs) const;
    */
  };
}

#endif
