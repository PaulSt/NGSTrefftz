#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include <l2hofefo.hpp>
#include "helpers.cpp"
// using namespace ngfem;

namespace ngfem
{
  template <int D, int ord> class TrefftzElement : public FiniteElement
  {
  private:
    constexpr static int nbasis
        = BinCoeff (D + ord, ord) + BinCoeff (D + ord - 1, ord - 1);

    constexpr static int npoly = BinCoeff (D + 1 + ord, ord);

    static const Mat<npoly, D + 1, int> indices;

    static const Mat<nbasis, npoly, double> basis;

  public:
    TrefftzElement () : FiniteElement ()
    {
      // cout << "ord: " + to_string(ord) + ", dimension: " + to_string(D) + ",
      // number of basis functions: " << nbasis << endl;
    }

    virtual ELEMENT_TYPE ElementType () const { return ET_TRIG; }

    virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                            Vector<double> &shape) const;

    virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                             SliceMatrix<> dshape) const;

    double ipow_ar (FlatVector<double> base, Vec<D + 1, int> ex,
                    float result = 1, int count = D + 1) const;

    int GetNBasis () const;

    void static MakeIndices_inner (Mat<npoly, D + 1, int> &indice,
                                   Vec<D + 1, int> &numbers, int &count,
                                   int dim = D + 1)
    {
      if (dim > 0)
        {
          for (int i = 0; i <= ord; i++)
            {
              numbers (D + 1 - dim) = i;
              MakeIndices_inner (indice, numbers, count, dim - 1);
            }
        }
      else
        {
          int sum = 0;
          for (int i = 0; i < D + 1; i++)
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

    constexpr static Mat<npoly, D + 1, int> MakeIndices ()
    {
      Mat<npoly, D + 1, int> indice = 0;
      Vec<D + 1, int> numbers = 0;
      int count = 0;
      MakeIndices_inner (indice, numbers, count);
      return indice;
    }

    constexpr static int IndexMap (Vec<D + 1, int> index)
    {
      int sum = 0;
      int temp_size = 0;
      for (int d = 0; d < D + 1; d++)
        {
          for (int p = 0; p < index (d); p++)
            {
              sum += BinCoeff (D - d + ord - p - temp_size,
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
               i++) // loop over indices BinCoeff(D+1 + ord, ord)
            {
              int k = indices (i, 0);
              if (k > 1)
                {
                  for (int m = 1; m <= D; m++) // rekursive sum
                    {
                      Vec<D + 1, int> get_coeff = indices.Row (i);
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

    // using ScalarFiniteElement<2>::CalcShape;
    // using ScalarFiniteElement<2>::CalcDShape;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement (py::module m);
#endif // NGS_PYTHON

#endif // FILE_TrefftzElement_HPP
