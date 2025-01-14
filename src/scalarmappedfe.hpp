#ifndef FILE_SCALARMAPPEDELEMENT_HPP
#define FILE_SCALARMAPPEDELEMENT_HPP

#include <fem.hpp>
#include "ngsttd.hpp"

namespace ngfem
{
  class Monomial : public RecursivePolynomial<Monomial>
  {
  public:
    Monomial () { ; }

    template <class S, class T> inline Monomial (int n, S x, T &&values)
    {
      Eval (n, x, values);
    }

    template <class S> static INLINE double P0 (S) { return 1.0; }
    template <class S> static INLINE S P1 (S x) { return x; }
    template <class S, class Sy> static INLINE S P1 (S x, Sy)
    {
      return P1 (x);
    }

    static INLINE double A (int) { return 1.0; }
    static INLINE double B (int) { return 0; }
    static INLINE double C (int) { return 0; }

    static INLINE double CalcA (int) { return 1.0; }
    static INLINE double CalcB (int) { return 0; }
    static INLINE double CalcC (int) { return 0; }

    enum
    {
      ZERO_B = 1
    };
  };

  template <int D, typename T> T vsum (Vec<D, T> v)
  {
    T sum = 0;
    for (int i = 0; i < D; i++)
      sum += v[i];
    return sum;
  }

  template <int D, typename T, typename T2>
  void vtimes (Vec<D, T> &v, Vec<D, T2> w)
  {
    for (int i = 0; i < D; i++)
      v[i] *= w[i];
  }

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
    NGST_DLL virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceVector<> shape) const
        = 0;

    NGST_DLL virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceVector<Complex> shape) const;

    // compute dshape, matrix: ndof x spacedim
    NGST_DLL virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                                      BareSliceMatrix<> dshape) const
        = 0;

    // returns shape functions in point ip.
    INLINE FlatVector<>
    GetShape (const BaseMappedIntegrationPoint &mip, LocalHeap &lh) const
    {
      FlatVector<> shape (ndof, lh);
      CalcShape (mip, shape);
      return shape;
    }

    // compute shape, row is shape nr, col is ip nr
    NGST_DLL virtual void CalcShape (const BaseMappedIntegrationRule &mir,
                                     BareSliceMatrix<> shape) const;

    // compute shape, row is shape nr, col is ip nr
    NGST_DLL virtual void
    CalcShape (const SIMD_BaseMappedIntegrationRule &mir,
               BareSliceMatrix<SIMD<double>> shape) const;

    NGST_DLL virtual void
    CalcMappedDShape (const SIMD_BaseMappedIntegrationRule &mir,
                      BareSliceMatrix<SIMD<double>> dshapes) const;

    // Evaluates function in integration point ip / integration rule ir.
    // Vector x provides coefficient vector.
    NGST_DLL virtual double Evaluate (const BaseMappedIntegrationPoint &mip,
                                      BareSliceVector<> x) const;
    NGST_DLL virtual void
    Evaluate (const BaseMappedIntegrationRule &mir, BareSliceVector<> coefs,
              FlatVector<> values) const;
    NGST_DLL virtual void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir,
              BareSliceVector<> coefs, BareVector<SIMD<double>> values) const;
    NGST_DLL virtual void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
              BareSliceMatrix<SIMD<double>> values) const;
    // Each column a vector ...
    NGST_DLL virtual void
    Evaluate (const BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
              SliceMatrix<> values) const;

    // Evaluate function in points of integrationrule ir, transpose operation.
    // Vector x provides coefficient vector.
    NGST_DLL virtual void EvaluateTrans (const BaseMappedIntegrationRule &mir,
                                         FlatVector<double> values,
                                         BareSliceVector<double> coefs) const;
    NGST_DLL virtual void
    AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
              BareVector<SIMD<double>> values, BareSliceVector<> coefs) const;
    NGST_DLL virtual void
    AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
              BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const;

    NGST_DLL virtual void
    EvaluateGrad (const SIMD_BaseMappedIntegrationRule &ir,
                  BareSliceVector<> coefs,
                  BareSliceMatrix<SIMD<double>> values) const;
    // needed for ALE-trafo
    // NGST_DLL virtual void EvaluateGrad (const SIMD_IntegrationRule
    // & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values)
    // const;
    NGST_DLL virtual void
    AddGradTrans (const SIMD_BaseMappedIntegrationRule &ir,
                  BareSliceMatrix<SIMD<double>> values,
                  BareSliceVector<> coefs) const;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef Vec<3, Array<double>>
      CSR; // CSR sparse matrix in (row,col,val) format

  void MatToCSR (Matrix<> mat, CSR &sparsemat);

  constexpr inline int BinCoeff (int n, int k) noexcept
  {
    return (k > n) ? 0 : // out of range
               (k == 0 || k == n) ? 1
                                  : // edge
               (k == 1 || k == n - 1) ? n
                                      : // first
               (k + k < n) ?            // recursive:
               (BinCoeff (n - 1, k - 1) * n) / k
                           :                        //  path to k=1   is faster
               (BinCoeff (n - 1, k) * n) / (n - k); //  path to k=n-1 is faster
  }

  constexpr long ipow (int base, int expo, int result = 1)
  {
    return expo < 1 ? result
                    : ipow (base * base, expo / 2,
                            (expo % 2) ? result * base : result);
  }

  template <int D> class ScalarMappedElement : public BaseScalarMappedElement
  {

  protected:
    CSR localmat;
    ELEMENT_TYPE eltype;
    Vec<D> shift;
    Vec<D> scale;
    const int npoly = BinCoeff (D + order, order);

  public:
    using BaseScalarMappedElement::BaseScalarMappedElement;

    ScalarMappedElement (int andof, int aord, CSR alocalmat,
                         ELEMENT_TYPE aeltype, Vec<D> ashift = 0,
                         Vec<D> ascale = 1.0)
        : BaseScalarMappedElement (andof, aord), localmat (alocalmat),
          eltype (aeltype), shift (ashift), scale (ascale)
    {
      ;
    }

    Vec<D> GetScale () const { return scale; }
    void SetScale (Vec<D> ascale) { scale = ascale; }

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
    // NGST_DLL virtual void GetPolOrders (FlatArray<PolOrder<D> >
    // orders) const;

    // public:
    //	NGST_DLL virtual std::list<std::tuple<std::string,double>> Timing
    //() const;

    // NGST_DLL virtual void CalcMappedDDShape (const
    // BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> hddshape) const;
    // using ScalarMappedElement<D+1>::CalcMappedDDShape;
    virtual void CalcMappedDDShape (const BaseMappedIntegrationPoint &bmip,
                                    BareSliceMatrix<> hddshape) const;
    void CalcDDWaveOperator (const SIMD_BaseMappedIntegrationRule &smir,
                             BareSliceMatrix<SIMD<double>> dshape,
                             BareSliceMatrix<SIMD<double>> wavespeed,
                             BareSliceMatrix<SIMD<double>> mu) const;

    void CalcDDWaveOperator (const SIMD_BaseMappedIntegrationRule &smir,
                             BareSliceMatrix<SIMD<double>> dshape,
                             BareSliceMatrix<SIMD<double>> wavespeed) const
    {
      Matrix<SIMD<double>> mu (1, wavespeed.Dist ());
      SIMD<double> a = 1.0;
      mu = a;
      CalcDDWaveOperator (smir, dshape, wavespeed, mu);
    }

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

  template <int D> class BlockMappedElement : public ScalarMappedElement<D>
  {
  private:
    Vector<CSR> localmats;

  public:
    BlockMappedElement (int andof, int aord, Vector<CSR> alocalmats,
                        ELEMENT_TYPE aeltype, Vec<D> ashift = 0,
                        Vec<D> ascale = 1.0)
        : ScalarMappedElement<D> (andof, aord, alocalmats[0], aeltype, ashift,
                                  ascale),
          localmats (alocalmats)
    {
      ;
    }

    // virtual ELEMENT_TYPE ElementType() const { return eltype; }

    // using ScalarMappedElement<D>::CalcShape;
    // using ScalarMappedElement<D>::CalcDShape;
    using ScalarMappedElement<D>::CalcMappedDShape;
    using ScalarMappedElement<D>::Evaluate;
    using ScalarMappedElement<D>::EvaluateGrad;
    using ScalarMappedElement<D>::AddGradTrans;

    void CalcShape (const BaseMappedIntegrationPoint &mip,
                    BareSliceVector<> shape) const override;
    void CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                    BareSliceMatrix<SIMD<double>> shape) const override;

    void CalcDShape (const BaseMappedIntegrationPoint &mip,
                     BareSliceMatrix<> dshape) const override;
    // void CalcDShape (const BaseMappedIntegrationRule & mir,
    // BareSliceMatrix<> dshapes) const override;
    void CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                     BareSliceMatrix<SIMD<double>> dshape) const override;
    void CalcMappedDDShape (const BaseMappedIntegrationPoint &bmip,
                            BareSliceMatrix<> hddshape) const override;
  };

}
#endif
