#ifndef FILE_SCALARMAPPEDELEMENT_HPP
#define FILE_SCALARMAPPEDELEMENT_HPP

#include <fem.hpp>

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
               BareSliceVector<> shape) const
        = 0;

    HD NGS_DLL_HEADER virtual void
    CalcShape (const BaseMappedIntegrationPoint &mip,
               BareSliceVector<Complex> shape) const;

    // compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER virtual void
    CalcDShape (const BaseMappedIntegrationPoint &mip,
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

  typedef Vec<3, Array<double>>
      CSR; // CSR sparse matrix in (row,col,val) format

  void MatToCSR (Matrix<> mat, CSR &sparsemat);

  constexpr inline size_t BinCoeff (size_t n, size_t k) noexcept
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
    Vec<D> elcenter;
    double elsize;
    float c = 1.0;
    const int npoly = BinCoeff (D + order, order);

  public:
    using BaseScalarMappedElement::BaseScalarMappedElement;

    ScalarMappedElement (int andof, int aord, CSR alocalmat,
                         ELEMENT_TYPE aeltype, Vec<D> aelcenter = 0,
                         double aelsize = 1, double ac = 1.0)
        : BaseScalarMappedElement (andof, aord), localmat (alocalmat),
          eltype (aeltype), elcenter (aelcenter), elsize (aelsize), c (ac)
    {
      ;
    }

    double GetWavespeed () const { return this->c; }
    void SetWavespeed (double wavespeed) { this->c = wavespeed; }

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
    // consistancy, can use CalcDShape with BaseMappedIR HD NGS_DLL_HEADER
    // virtual void CalcMappedDShape (const BaseMappedIntegrationPoint & mip,
    // BareSliceMatrix<> dshape) const;
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

    // NGS_DLL_HEADER virtual void CalcMappedDDShape (const
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
    // const int nbasis;
    // const int npoly;
    // Vec<D> elcenter;
    // double elsize;
    // float c;
    // ELEMENT_TYPE eltype;
    Vector<CSR> localmats;

  public:
    // TrefftzWaveFE();
    BlockMappedElement (int andof, int aord, Vector<CSR> alocalmats,
                        ELEMENT_TYPE aeltype, Vec<D> aelcenter = 0,
                        double aelsize = 1, float ac = 1.0)
        : ScalarMappedElement<D> (andof, aord, alocalmats[0], aeltype,
                                  aelcenter, aelsize, ac),
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
