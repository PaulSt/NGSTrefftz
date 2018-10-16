#include "trefftzelement.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
  template <int D>
  T_TrefftzElement<D>::T_TrefftzElement (int aord, float ac,
                                         ELEMENT_TYPE aeltype, int abasistype)
      : ScalarMappedElement<D> (BinCoeff (D - 1 + aord, aord)
                                    + BinCoeff (D - 1 + aord - 1, aord - 1),
                                aord),
        ord (aord), c (ac), nbasis (BinCoeff (D - 1 + ord, ord)
                                    + BinCoeff (D - 1 + ord - 1, ord - 1)),
        npoly (BinCoeff (D + ord, ord)), indices (GetIndices ()),
        basistype (abasistype), eltype (aeltype), pascal (pascal_sym ())
  {
    ;
  }

  template <int D>
  void T_TrefftzElement<D>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                       BareSliceVector<> shape) const
  {
    Vec<D> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[D - 1] *= c;
    Matrix<double> coeff (TrefftzBasis ());

    for (int j = ord; j > 0; j--)
      for (int i = 0; i < D; i++)
        for (int k = pascal (i + 1, j) - 1; k >= 0; k--)
          coeff.Row (pascal (D + 1, j - 1) + k)
              += cpoint[i]
                 * coeff.Row (pascal (D + 1, j) + pascal (i, j + 1) + k);

    for (int b = 0; b < nbasis; b++)
      shape (b) = coeff.Row (0) (b);
  }

  template <int D>
  void T_TrefftzElement<D>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                        SliceMatrix<> dshape) const
  {
    Vec<D> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[D - 1] *= c;

    for (int d = 0; d < D; d++)
      { // loop over derivatives/dimensions
        Matrix<double> coeff (GetDerTrefftzBasis (d));
        for (int j = ord - 1; j > 0; j--)
          for (int i = 0; i < D; i++)
            for (int k = pascal (i + 1, j) - 1; k >= 0; k--)
              coeff.Row (pascal (D + 1, j - 1) + k)
                  += cpoint[i]
                     * coeff.Row (pascal (D + 1, j) + pascal (i, j + 1) + k);
        dshape.Col (d) = coeff.Row (0);
      }

    dshape.Col (D - 1) *= c;  // inner derivative
    dshape *= (2.0 / elsize); // inner derivative
  }

  template <int D>
  Matrix<double> T_TrefftzElement<D>::GetDerTrefftzBasis (int der) const
  {
    static int order;
    static int btype;
    static Vec<D, Matrix<double>> basisstorage;
    if (order != ord || btype != basistype)
      {
        for (int d = 0; d < D; d++)
          {
            basisstorage[d].SetSize (pascal (D + 1, ord), nbasis);
            for (int i = 0, count = 0; i < npoly; i++)
              {
                if (indices (i, d) != 0)
                  basisstorage[d].Row (count++)
                      = indices (i, d) * TrefftzBasis ().Row (i);
              }
          }
        order = ord;
        btype = basistype;
      }
    return basisstorage[der];
  }

  template <int D> Matrix<double> T_TrefftzElement<D>::TrefftzBasis () const
  {
    static int order;
    static int btype;
    static Matrix<double> basisstorage;

    if (order != ord || btype != basistype)
      {
        basisstorage.SetSize (npoly, nbasis);

        basisstorage = 0;
        int setbasis = 0;
        for (int l = 0; l < nbasis; l++) // loop over basis functions
          {
            for (int i = 0; i < npoly;
                 i++) // loop over indices BinCoeff(D + ord, ord)
              {
                int k = indices (i, D - 1);
                if (k > 1)
                  {
                    for (int m = 0; m < D - 1; m++) // rekursive sum
                      {
                        Vec<D, int> get_coeff = indices.Row (i);
                        get_coeff[D - 1] = get_coeff[D - 1] - 2;
                        get_coeff[m] = get_coeff[m] + 2;
                        basisstorage (i, l)
                            += (indices (i, m) + 1) * (indices (i, m) + 2)
                               * basisstorage (IndexMap (get_coeff), l);
                      }
                    basisstorage (i, l) *= 1.0 / (k * (k - 1));
                  }
                else if (k <= 1) // time=0 and =1
                  {
                    switch (basistype)
                      {
                      case 0:
                        if (l == 0)
                          basisstorage (i, setbasis++)
                              = 1.0; // set the l-th coeff to 1
                        // i += nbasis-1;	//jump to time = 2 if i=0
                        break;
                      case 1:
                        if ((k == 0 && l < BinCoeff (D - 1 + ord, ord))
                            || (k == 1 && l >= BinCoeff (D - 1 + ord, ord)))
                          {
                            basisstorage (i, l) = 1;
                            for (int exponent :
                                 indices.Row (i).Range (0, D - 1))
                              basisstorage (i, l)
                                  *= LegCoeffMonBasis (l, exponent);
                          }
                        break;
                      case 2:
                        if ((k == 0 && l < BinCoeff (D - 1 + ord, ord))
                            || (k == 1 && l >= BinCoeff (D - 1 + ord, ord)))
                          {
                            basisstorage (i, l) = 1;
                            for (int exponent :
                                 indices.Row (i).Range (0, D - 1))
                              basisstorage (i, l)
                                  *= ChebCoeffMonBasis (l, exponent);
                          }
                        break;
                      }
                  }
              }
          }
        order = ord;
        btype = basistype;
      }
    return basisstorage;
  }

  template <int D>
  constexpr void
  T_TrefftzElement<D>::MakeIndices_inner (Matrix<int> &indice,
                                          Vec<D, int> numbers, int &count,
                                          int ordr, int dim)
  {
    if (dim > 0)
      {
        for (int i = 0; i <= ordr; i++)
          {
            numbers (dim - 1) = i;
            MakeIndices_inner (indice, numbers, count, ordr, dim - 1);
          }
      }
    else
      {
        int sum = 0;
        for (int i = 0; i < D; i++)
          {
            sum += numbers (i);
          }
        if (sum == ordr)
          {
            indice.Row (count++) = numbers;
          }
      }
  }

  template <int D> Matrix<int> T_TrefftzElement<D>::GetIndices ()
  {
    static int order;
    static Matrix<int> indice;
    if (order != ord)
      {
        indice.SetSize (npoly, D);
        Vec<D, int> numbers = 0;
        int count = 0;
        for (int o = 0; o <= ord; o++)
          {
            MakeIndices_inner (indice, numbers, count, o);
          }
      }
    return indice;
  }

  template <int D>
  constexpr int T_TrefftzElement<D>::IndexMap (Vec<D, int> index) const
  {
    int sum = 0;
    int indexleng = 0;
    for (int r = 0; r < D; r++)
      {
        indexleng += index (r);
        for (int i = 0; i < index (r); i++)
          {
            sum += BinCoeff (indexleng - i + r - 1, indexleng - i);
          }
      }
    sum += BinCoeff (indexleng - 1 + D, indexleng - 1);
    return sum;
  }

  template <int D> Matrix<int> T_TrefftzElement<D>::pascal_sym () const
  {
    static int order;
    static Matrix<int> pascalstorage;

    if (order != ord)
      {
        pascalstorage.SetSize (D + 2, ord + 2);
        for (int i = 0; i <= D + 1; ++i)
          for (int j = 0; j <= ord + 1; ++j)
            if (i == 0 || j == 0)
              pascalstorage (i, j) = 0;
            else if (i == 1 || j == 1)
              pascalstorage (i, j) = 1;
            else
              pascalstorage (i, j)
                  = pascalstorage (i - 1, j) + pascalstorage (i, j - 1);

        order = ord;
      }
    return pascalstorage;
  }

  template class T_TrefftzElement<1>;
  template class T_TrefftzElement<2>;
  template class T_TrefftzElement<3>;
  template class T_TrefftzElement<4>;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NGS_PYTHON
void ExportTrefftzElement (py::module m)
{
  // py::class_<T_TrefftzElement<3>, shared_ptr<T_TrefftzElement<3>>,
  // FiniteElement> 	(m, "T_TrefftzElement3", "Trefftz space for wave eq")
  // 	.def(py::init<>())
  // 	;
  // py::class_<T_TrefftzElement<2>, shared_ptr<T_TrefftzElement<2>>,
  // FiniteElement> 	(m, "T_TrefftzElement2", "Trefftz space for wave eq")
  // 	.def(py::init<>())
  // 	;
  // py::class_<T_TrefftzElement<1>, shared_ptr<T_TrefftzElement<1>>,
  // FiniteElement> 	(m, "T_TrefftzElement1", "Trefftz space for wave eq")
  // 	.def(py::init<>())
  // 	;
}
#endif // NGS_PYTHON
