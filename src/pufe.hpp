#ifndef FILE_PUFELEMENT_HPP
#define FILE_PUFELEMENT_HPP

#include <fem.hpp>
#include "ngsttd.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{

  template <int D> class PUFElement : public BaseScalarMappedElement
  {

  protected:
    CSR localmat;
    ELEMENT_TYPE eltype;
    Vec<D + 1, Vec<D>> elvertices;
    Vec<D + 1> elsizes;
    float c = 1.0;
    const int npoly = BinCoeff (D + order, order);

  public:
    using BaseScalarMappedElement::BaseScalarMappedElement;

    PUFElement (int andof, int aord, CSR alocalmat, ELEMENT_TYPE aeltype,
                Vec<D + 1, Vec<D>> aelvertices = 0, Vec<D + 1> aelsizes = 1,
                double ac = 1.0)
        : BaseScalarMappedElement (andof, aord), localmat (alocalmat),
          eltype (aeltype), elvertices (aelvertices), elsizes (aelsizes),
          c (ac)
    {
      ;
    }

    double GetWavespeed () const { return this->c; }
    void SetWavespeed (double wavespeed) { this->c = wavespeed; }

    // the name
    NGST_DLL virtual string ClassName () const;
    NGST_DLL virtual int Dim () const { return D; }

    // returns derivatives in point ip.
    INLINE const FlatMatrixFixWidth<D>
    GetDShape (const BaseMappedIntegrationPoint &mip, LocalHeap &lh) const
    {
      FlatMatrixFixWidth<D> dshape (ndof, lh);
      CalcDShape (mip, dshape);
      return dshape;
    }

    virtual ELEMENT_TYPE ElementType () const { return eltype; }

    using BaseScalarMappedElement::CalcDShape;
    using BaseScalarMappedElement::CalcMappedDShape;
    using BaseScalarMappedElement::CalcShape;

    virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                            BareSliceVector<> shape) const;
    virtual void CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                            BareSliceMatrix<SIMD<double>> shape) const;
    virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                             BareSliceMatrix<> dshape) const;
    virtual void CalcDShape (const BaseMappedIntegrationRule &mir,
                             BareSliceMatrix<> dshapes) const;
    virtual void CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                             BareSliceMatrix<SIMD<double>> dshape) const;

    // compute dshape, matrix: ndof x spacedim, Use CalcMappedDShape only for
    // consistancy, can use CalcDShape with BaseMappedIR NGST_DLL
    // virtual void CalcMappedDShape (const BaseMappedIntegrationPoint & mip,
    // BareSliceMatrix<> dshape) const;
    NGST_DLL virtual void
    CalcMappedDShape (const BaseMappedIntegrationRule &mir,
                      BareSliceMatrix<> dshapes) const;

    // Evaluates gradient in integration point ip.
    // Vector x provides coefficient vector.
    NGST_DLL virtual Vec<D> EvaluateGrad (const BaseMappedIntegrationPoint &ip,
                                          BareSliceVector<> x) const;

    using BaseScalarMappedElement::AddGradTrans;
    using BaseScalarMappedElement::Evaluate;
    using BaseScalarMappedElement::EvaluateGrad;

    // Evaluate gradient in points of integrationrule ir.
    // Vector x provides coefficient vector.
    NGST_DLL virtual void
    EvaluateGrad (const BaseMappedIntegrationRule &ir, BareSliceVector<> coefs,
                  FlatMatrixFixWidth<D> values) const;

    // Evaluate gradient in points of integrationrule ir, transpose operation.
    // Vector x provides coefficient vector.
    NGST_DLL virtual void
    EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                       FlatMatrixFixWidth<D> values,
                       BareSliceVector<> coefs) const;
    NGST_DLL virtual void
    EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                       SliceMatrix<> values, SliceMatrix<> coefs) const;

    void CalcMappedDShape (const BaseMappedIntegrationPoint &bmip,
                           BareSliceMatrix<> dshape) const
    {
      CalcDShape (bmip, dshape);
    }
    void CalcMappedDShape (const SIMD_BaseMappedIntegrationRule &mir,
                           BareSliceMatrix<SIMD<double>> dshapes) const
    {
      CalcDShape (mir, dshapes);
    }

    void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir,
              BareSliceVector<> coefs, BareVector<SIMD<double>> values) const
    {
      STACK_ARRAY (SIMD<double>, mem, this->ndof * mir.Size ());
      FlatMatrix<SIMD<double>> shape (this->ndof, mir.Size (), &mem[0]);
      CalcShape (mir, shape);
      const int nsimd = SIMD<double>::Size ();
      FlatMatrix<double> bdbmat (this->ndof, mir.Size () * nsimd,
                                 reinterpret_cast<double *> (&shape (0, 0)));
      FlatVector<double> bdbvec (mir.Size () * nsimd,
                                 reinterpret_cast<double *> (&values (0)));
      bdbvec = Trans (bdbmat) * coefs;
    }
    void
    AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
              BareVector<SIMD<double>> values, BareSliceVector<> coefs) const
    {
      STACK_ARRAY (SIMD<double>, mem, this->ndof * mir.Size ());
      FlatMatrix<SIMD<double>> shape (this->ndof, mir.Size (), &mem[0]);
      CalcShape (mir, shape);
      const int nsimd = SIMD<double>::Size ();
      FlatMatrix<double> bdbmat (this->ndof, mir.Size () * nsimd,
                                 reinterpret_cast<double *> (&shape (0, 0)));
      FlatVector<double> bdbvec (mir.Size () * nsimd,
                                 reinterpret_cast<double *> (&values (0)));
      coefs.Range (0, this->ndof) += bdbmat * bdbvec;
    }

    void EvaluateGrad (const SIMD_BaseMappedIntegrationRule &ir,
                       BareSliceVector<> coefs,
                       BareSliceMatrix<SIMD<double>> values) const
    {
      STACK_ARRAY (SIMD<double>, mem, D * this->ndof * ir.Size ());
      FlatMatrix<SIMD<double>> simddshapes (D * this->ndof, ir.Size (),
                                            &mem[0]);
      CalcDShape (ir, simddshapes);
      const int nsimd = SIMD<double>::Size ();
      FlatMatrix<double> dshapes (
          this->ndof, D * nsimd * ir.Size (),
          reinterpret_cast<double *> (&simddshapes (0, 0)));
      FlatVector<double> bdbvec (D * nsimd * ir.Size (),
                                 reinterpret_cast<double *> (&values (0, 0)));
      bdbvec = Trans (dshapes) * coefs;
    }
    void AddGradTrans (const SIMD_BaseMappedIntegrationRule &mir,
                       BareSliceMatrix<SIMD<double>> values,
                       BareSliceVector<> coefs) const
    {
      STACK_ARRAY (SIMD<double>, mem, D * this->ndof * mir.Size ());
      FlatMatrix<SIMD<double>> simddshapes (D * this->ndof, mir.Size (),
                                            &mem[0]);
      CalcDShape (mir, simddshapes);
      const int nsimd = SIMD<double>::Size ();
      FlatMatrix<double> dshapes (
          this->ndof, D * nsimd * mir.Size (),
          reinterpret_cast<double *> (&simddshapes (0, 0)));
      FlatVector<double> bdbvec (D * nsimd * mir.Size (),
                                 reinterpret_cast<double *> (&values (0, 0)));
      coefs.Range (0, this->ndof) += dshapes * bdbvec;
    }
  };

}
#endif
