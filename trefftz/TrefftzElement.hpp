#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include <l2hofefo.hpp>
#include "helpers.cpp"
#include "MappedElement.hpp"
// using namespace ngfem;

namespace ngfem
{
  template <int D, int ord> class TrefftzElement : public MappedElement
  {
  private:
    constexpr static int nbasis
        = BinCoeff (D - 1 + ord, ord) + BinCoeff (D - 1 + ord - 1, ord - 1);

    constexpr static int npoly = BinCoeff (D + ord, ord);

    static const Mat<npoly, D, int> indices;

    // static const Mat<nbasis, npoly,double> basis;
    Matrix<double> basis;

  public:
    TrefftzElement () : MappedElement (), basis (npoly, D)
    {
      basis = TrefftzBasis ();
      // cout << "ord: " + to_string(ord) + ", dimension: " + to_string(D) + ",
      // number of basis functions: " << nbasis << endl;
    }

    virtual ELEMENT_TYPE ElementType () const { return ET_TRIG; }

    virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                            BareSliceVector<> shape) const;

    virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                             SliceMatrix<> dshape) const;

    double ipow_ar (FlatVector<double> base, Vec<D, int> ex, float result = 1,
                    int count = D) const;

    int GetNBasis () const;

    void static MakeIndices_inner (Mat<npoly, D, int> &indice,
                                   Vec<D, int> &numbers, int &count,
                                   int dim = D)
    {
      if (dim > 0)
        {
          for (int i = 0; i <= ord; i++)
            {
              numbers (D - dim) = i;
              MakeIndices_inner (indice, numbers, count, dim - 1);
            }
        }
      else
        {
          int sum = 0;
          for (int i = 0; i < D; i++)
            {
              sum += numbers (i);
            }
          if (sum <= ord)
            {
              indice.Row (count++) = numbers;
              // cout << IndexMap(indices.Row(count-1)) << ": " <<
              // indices.Row(count-1) << endl;
            }
        }
    }

    constexpr static Mat<npoly, D, int> MakeIndices ()
    {
      Mat<npoly, D, int> indice = 0;
      Vec<D, int> numbers = 0;
      int count = 0;
      MakeIndices_inner (indice, numbers, count);
      return indice;
    }

    constexpr static int IndexMap (Vec<D, int> index)
    {
      int sum = 0;
      int temp_size = 0;
      for (int d = 0; d < D; d++)
        {
          for (int p = 0; p < index (d); p++)
            {
              sum += BinCoeff (D - 1 - d + ord - p - temp_size,
                               ord - p - temp_size);
            }
          temp_size += index (d);
        }
      return sum;
    }

    constexpr static Mat<nbasis, npoly, double> TrefftzBasis ()
    {
      Mat<nbasis, npoly, double> temp_basis = 0;
      for (int l = 0; l < nbasis; l++) // loop over basis functions
        {
          for (int i = 0; i < npoly;
               i++) // loop over indices BinCoeff(D + ord, ord)
            {
              int k = indices (i, 0);
              if (k > 1)
                {
                  for (int m = 1; m <= D - 1; m++) // rekursive sum
                    {
                      Vec<D, int> get_coeff = indices.Row (i);
                      get_coeff[0] = get_coeff[0] - 2;
                      get_coeff[m] = get_coeff[m] + 2;
                      temp_basis (l, IndexMap (indices.Row (i)))
                          += (indices (i, m) + 1) * (indices (i, m) + 2)
                             * temp_basis (l, IndexMap (get_coeff));
                    }
                  temp_basis (l, IndexMap (indices.Row (i)))
                      *= 1.0 / (k * (k - 1));
                }
              else if (k == 0) // time=0
                {
                  temp_basis (l, IndexMap (indices.Row (l))) = 1.0;
                  i += nbasis;
                }
            }
        }
      return temp_basis;
      // cout << "basis: \n" << basis << endl;
    }

    using MappedElement::CalcDShape;
    using MappedElement::CalcShape;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement (py::module m);
#endif // NGS_PYTHON

#endif // FILE_TrefftzElement_HPP
