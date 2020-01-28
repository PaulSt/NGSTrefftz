#ifndef FILE_SCALARMAPPEDELEMENT_HPP
#define FILE_SCALARMAPPEDELEMENT_HPP

#include <fem.hpp>
#include "helpers.hpp"

namespace ngfem
{

  class BaseScalarMappedElement : public FiniteElement
  {
  public:
    // using FiniteElement::FiniteElement;

    INLINE BaseScalarMappedElement () { ; }
    INLINE BaseScalarMappedElement (int andof, int aorder)
        : FiniteElement (andof, aorder)
    {
      ;
    }

    // compute shape
    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationPoint &mip,
               BareSliceVector<> shape) const = 0;

    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationPoint &mip,
               BareSliceVector<Complex> shape) const;

    // compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER virtual void
    CalcDShape (const BaseMappedIntegrationPoint &mip,
                BareSliceMatrix<> dshape) const = 0;

    // returns shape functions in point ip.
    INLINE FlatVector<>
    GetShape (const BaseMappedIntegrationPoint &mip, LocalHeap &lh) const
    {
      FlatVector<> shape (ndof, lh);
      CalcShape (mip, shape);
      return shape;
    }

    // compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationRule &mir,
               SliceMatrix<> shape) const;

    // compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER virtual void
    CalcShape (const SIMD_BaseMappedIntegrationRule &mir,
               BareSliceMatrix<SIMD<double>> shape) const;

    HD NGS_DLL_HEADER virtual void
    CalcMappedDShape (const SIMD_BaseMappedIntegrationRule &mir,
                      BareSliceMatrix<SIMD<double>> dshapes) const;

    // Evaluates function in integration point ip / integration rule ir.
    // Vector x provides coefficient vector.
    HD NGS_DLL_HEADER virtual double
    Evaluate (const BaseMappedIntegrationPoint &mip,
              BareSliceVector<> x) const;
    HD NGS_DLL_HEADER virtual void
    Evaluate (const BaseMappedIntegrationRule &mir, BareSliceVector<> coefs,
              FlatVector<> values) const;
    HD NGS_DLL_HEADER virtual void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir,
              BareSliceVector<> coefs, BareVector<SIMD<double>> values) const;
    HD NGS_DLL_HEADER virtual void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
              BareSliceMatrix<SIMD<double>> values) const;
    // Each column a vector ...
    HD NGS_DLL_HEADER virtual void
    Evaluate (const BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
              SliceMatrix<> values) const;

    // Evaluate function in points of integrationrule ir, transpose operation.
    // Vector x provides coefficient vector.
    HD NGS_DLL_HEADER virtual void
    EvaluateTrans (const BaseMappedIntegrationRule &mir,
                   FlatVector<double> values,
                   BareSliceVector<double> coefs) const;
    HD NGS_DLL_HEADER virtual void
    AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
              BareVector<SIMD<double>> values, BareSliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void
    AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
              BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const;

    HD NGS_DLL_HEADER virtual void
    EvaluateGrad (const SIMD_BaseMappedIntegrationRule &ir,
                  BareSliceVector<> coefs,
                  BareSliceMatrix<SIMD<double>> values) const;
    // needed for ALE-trafo
    // HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_IntegrationRule
    // & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values)
    // const;
    HD NGS_DLL_HEADER virtual void
    AddGradTrans (const SIMD_BaseMappedIntegrationRule &ir,
                  BareSliceMatrix<SIMD<double>> values,
                  BareSliceVector<> coefs) const;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D> class ScalarMappedElement : public BaseScalarMappedElement
  {
  public:
    using BaseScalarMappedElement::BaseScalarMappedElement;

    // the name
    NGS_DLL_HEADER virtual string ClassName () const;
    HD NGS_DLL_HEADER virtual int Dim () const { return D; }

    // returns derivatives in point ip.
    INLINE const FlatMatrixFixWidth<D>
    GetDShape (const BaseMappedIntegrationPoint &mip, LocalHeap &lh) const
    {
      FlatMatrixFixWidth<D> dshape (ndof, lh);
      CalcDShape (mip, dshape);
      return dshape;
    }

    using BaseScalarMappedElement::CalcDShape;
    using BaseScalarMappedElement::CalcMappedDShape;
    using BaseScalarMappedElement::CalcShape;

    void CalcDShape (const BaseMappedIntegrationRule &mir,
                     BareSliceMatrix<> dshapes) const;
    virtual void CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                             BareSliceMatrix<SIMD<double>> dshape) const;

    // compute dshape, matrix: ndof x spacedim, Use CalcMappedDShape only for
    // consistancy, can use CalcDShape with BaseMappedIR
    HD NGS_DLL_HEADER virtual void
    CalcMappedDShape (const BaseMappedIntegrationPoint &mip,
                      BareSliceMatrix<> dshape) const;
    HD NGS_DLL_HEADER virtual void
    CalcMappedDShape (const BaseMappedIntegrationRule &mir,
                      SliceMatrix<> dshapes) const;

    // Evaluates gradient in integration point ip.
    // Vector x provides coefficient vector.
    HD NGS_DLL_HEADER virtual Vec<D>
    EvaluateGrad (const BaseMappedIntegrationPoint &ip,
                  BareSliceVector<> x) const;

    using BaseScalarMappedElement::AddGradTrans;
    using BaseScalarMappedElement::Evaluate;
    using BaseScalarMappedElement::EvaluateGrad;

    // Evaluate gradient in points of integrationrule ir.
    // Vector x provides coefficient vector.
    HD NGS_DLL_HEADER virtual void
    EvaluateGrad (const BaseMappedIntegrationRule &ir, BareSliceVector<> coefs,
                  FlatMatrixFixWidth<D> values) const;

    // Evaluate gradient in points of integrationrule ir, transpose operation.
    // Vector x provides coefficient vector.
    HD NGS_DLL_HEADER virtual void
    EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                       FlatMatrixFixWidth<D> values,
                       BareSliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void
    EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                       SliceMatrix<> values, SliceMatrix<> coefs) const;
    // HD NGS_DLL_HEADER virtual void GetPolOrders (FlatArray<PolOrder<D> >
    // orders) const;

    // public:
    //	NGS_DLL_HEADER virtual std::list<std::tuple<std::string,double>> Timing
    //() const;
    virtual float GetWavespeed () const
    {
      return 0;
    } // ugly parent hack for trefftzwave

    NGS_DLL_HEADER virtual void
    CalcMappedDDShape (const BaseMappedIntegrationPoint &bmip,
                       BareSliceMatrix<> hddshape) const;
  };

}
#endif
