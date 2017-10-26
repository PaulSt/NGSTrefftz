#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include "helpers.cpp"

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
    HD NGS_DLL_HEADER virtual void
    CalcMappedDShape (const SIMD_BaseMappedIntegrationRule &mir,
                      BareSliceMatrix<SIMD<double>> dshapes) const;

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
    HD NGS_DLL_HEADER virtual void
    Evaluate (const BaseMappedIntegrationRule &mir, BareSliceVector<> coefs,
              FlatVector<> values) const;
    HD NGS_DLL_HEADER virtual void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir,
              BareSliceVector<> coefs, BareVector<SIMD<double>> values) const;
    HD NGS_DLL_HEADER virtual void
    Evaluate (const SIMD_BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
              BareSliceMatrix<SIMD<double>> values) const;

    /**
       Each column a vector ...
     */
    HD NGS_DLL_HEADER virtual void
    Evaluate (const BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
              SliceMatrix<> values) const;

    /**
       Evaluate function in points of integrationrule ir, transpose operation.
       Vector x provides coefficient vector.
     */

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
    HD NGS_DLL_HEADER virtual void
    EvaluateGrad (const SIMD_IntegrationRule &ir, BareSliceVector<> coefs,
                  BareSliceMatrix<SIMD<double>> values) const;
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

    /// the name
    NGS_DLL_HEADER virtual string ClassName () const;

    HD NGS_DLL_HEADER virtual int Dim () const { return D; }

    /**
             returns derivatives in point ip.
    */

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

    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER virtual void
    CalcMappedDShape (const MappedIntegrationPoint<D, D> &mip,
                      SliceMatrix<> dshape) const;

    HD NGS_DLL_HEADER virtual void
    CalcMappedDShape (const MappedIntegrationRule<D, D> &mir,
                      SliceMatrix<> dshapes) const;

    /*
                             returns second derivatives in point ip.
                             returns stored values for valid ip.IPNr(), else
       computes values
                    /
                    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip,
       LocalHeap & lh) const
                    {
                            FlatMatrix<> ddshape(ndof, D*D, lh);
                            CalcDDShape (ip, ddshape);
                            return ddshape;
                    }

                    /// compute dshape, matrix: ndof x (spacedim spacedim)
                    NGS_DLL_HEADER virtual void CalcDDShape (const
       IntegrationPoint & ip, FlatMatrix<> ddshape) const;

                    /// compute dshape, matrix: ndof x (spacedim spacedim)
                    NGS_DLL_HEADER virtual void CalcMappedDDShape (const
       MappedIntegrationPoint<D,D> & mip, SliceMatrix<> ddshape) const;
    */

    /**
             Evaluates gradient in integration point ip.
             Vector x provides coefficient vector.
     */

    HD NGS_DLL_HEADER virtual Vec<D>
    EvaluateGrad (const BaseMappedIntegrationPoint &ip,
                  BareSliceVector<> x) const;

    using BaseScalarMappedElement::AddGradTrans;
    using BaseScalarMappedElement::Evaluate;
    using BaseScalarMappedElement::EvaluateGrad;

    /**
             Evaluate gradient in points of integrationrule ir.
             Vector x provides coefficient vector.
     */

    HD NGS_DLL_HEADER virtual void
    EvaluateGrad (const BaseMappedIntegrationRule &ir, BareSliceVector<> coefs,
                  FlatMatrixFixWidth<D> values) const;

    /**
             Evaluate gradient in points of integrationrule ir, transpose
       operation. Vector x provides coefficient vector.
     */

    HD NGS_DLL_HEADER virtual void
    EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                       FlatMatrixFixWidth<D> values,
                       BareSliceVector<> coefs) const;

    HD NGS_DLL_HEADER virtual void
    EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                       SliceMatrix<> values, SliceMatrix<> coefs) const;

    HD NGS_DLL_HEADER virtual void
    GetPolOrders (FlatArray<PolOrder<D>> orders) const;

    // public:
    //	NGS_DLL_HEADER virtual std::list<std::tuple<std::string,double>> Timing
    //() const;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D, int ord>
  class TrefftzElement : public ScalarMappedElement<D>
  {
  private:
    constexpr static int nbasis
        = BinCoeff (D - 1 + ord, ord) + BinCoeff (D - 1 + ord - 1, ord - 1);
    constexpr static int npoly = BinCoeff (D + ord, ord);
    static const Mat<npoly, D, int> indices;
    // static const Mat<nbasis, npoly,double> basis;
    static const Matrix<double> basis;
    Vec<D> elcenter = 0;
    float elsize = 1;
    float c = 1;

  protected:
    void static MakeIndices_inner (Mat<npoly, D, int> &indice,
                                   Vec<D, int> &numbers, int &count,
                                   int dim = D);
    constexpr static Mat<npoly, D, int> MakeIndices ();
    constexpr static int IndexMap (Vec<D, int> index);

  public:
    TrefftzElement () : ScalarMappedElement<D> (nbasis, ord)
    {
      ;
    } // BaseScalarMappedElement(nbasis,ord) { ;	}//

    virtual ELEMENT_TYPE ElementType () const { return ET_TRIG; }

    using ScalarMappedElement<D>::CalcShape;
    virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                            BareSliceVector<> shape) const;

    using ScalarMappedElement<D>::CalcDShape;
    virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                             SliceMatrix<> dshape) const;

    constexpr static Mat<nbasis, npoly, double> TrefftzBasis ();

    int GetNBasis () const { return nbasis; }

    TrefftzElement<D, ord> *SetCenter (Vec<D> acenter)
    {
      elcenter = acenter;
      return this;
    }
    TrefftzElement<D, ord> *SetElSize (float aelsize)
    {
      elsize = aelsize;
      return this;
    }
    TrefftzElement<D, ord> *SetWavespeed (float ac)
    {
      c = ac;
      return this;
    }

    FlatVector<double> ShiftPoint (FlatVector<double> point) const
    {
      point -= elcenter;
      point *= (1.0 / elsize);
      point (0) *= c;
      return point;
    }

    double ipow_ar (FlatVector<double> base, Vec<D, int> ex, float result = 1,
                    int count = D - 1) const;

    double ipowD_ar (int der, FlatVector<double> base, Vec<D, int> ex,
                     float result = 1, int count = D - 1) const;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement (py::module m);
#endif // NGS_PYTHON

#endif // FILE_TrefftzElement_HPP
