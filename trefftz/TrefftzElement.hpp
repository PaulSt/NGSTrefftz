#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include "MappedElement.hpp"
#include "helpers.cpp"

namespace ngfem
{
  template <int D, int ord>
  class TrefftzElement : public ScalarMappedElement<D>
  {
  private:
    constexpr static int nbasis
        = BinCoeff (D - 1 + ord, ord) + BinCoeff (D - 1 + ord - 1, ord - 1);

    constexpr static int npoly = BinCoeff (D + ord, ord);

    static const Mat<npoly, D, int> indices;

    Vec<D> elcenter;

    // static const Mat<nbasis, npoly,double> basis;
    static const Matrix<double> basis;

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

    using BaseScalarMappedElement::CalcShape;
    virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                            BareSliceVector<> shape) const;

    using BaseScalarMappedElement::CalcDShape;
    virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                             SliceMatrix<> dshape) const;

    constexpr static Mat<nbasis, npoly, double> TrefftzBasis ();

    int GetNBasis () const { return nbasis; }

    TrefftzElement<D, ord> *SetCenter (Vec<D> acenter)
    {
      elcenter = acenter;
      return this;
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
