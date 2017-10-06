#include <fem.hpp>
#include "TrefftzElement.hpp"
//#include "helpers.cpp"

namespace ngfem
{

  template <int D, int ord>
  const Mat<TrefftzElement<D, ord>::npoly, D, int>
      TrefftzElement<D, ord>::indices = MakeIndices ();

  // template<int D, int ord>
  // const Mat<TrefftzElement<D,ord>::nbasis,
  // TrefftzElement<D,ord>::npoly,double> TrefftzElement<D,ord> :: basis =
  // TrefftzBasis();

  template <int D, int ord>
  const Matrix<double> TrefftzElement<D, ord>::basis = TrefftzBasis ();

  template <int D, int ord>
  void
  TrefftzElement<D, ord>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceVector<> shape) const
  {
    Vec<npoly, float> polynomial;

    FlatVector<double> point = mip.GetPoint ();
    // point = point - elcenter;
    for (int i = 0; i < npoly; i++) // loop over indices
      {
        polynomial (i) = ipow_ar (point, indices.Row (i));
      }

    // Vec<nbasis,double> tempshape;
    Vector<double> tempshape (nbasis);
    tempshape = basis * polynomial;
    for (int b = 0; b < nbasis; b++)
      shape (b) = tempshape (b); // loop over basis TODO replace this by
                                 // correct way of filling  BareSliceVector
    // FlatVec<nbasis> (&shape(0)) = tempshape;
    // FlatVector<> (nbasis,&shape(0)) = tempshape;
  }

  template <int D, int ord>
  void
  TrefftzElement<D, ord>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                      SliceMatrix<> dshape) const
  {
    FlatVector<double> point = mip.GetPoint ();
    Vec<npoly, float> polynomial;
    Vector<double> tempshape (nbasis);
    Mat<npoly, D, int> derindices;
    for (int d = 0; d < D; d++) // loop over derivatives/dimensions
      {
        derindices = MakeIndices ();
        for (int i = 0; i < npoly; i++) // loop over indices
          {
            derindices (i, d) = derindices (i, d) - 1;
            polynomial (i) = ipowD_ar (d, point, derindices.Row (i));
          }
        dshape.Col (d) = basis * polynomial;
      }
  }

  template <int D, int ord>
  void TrefftzElement<D, ord>::MakeIndices_inner (
      Mat<TrefftzElement<D, ord>::npoly, D, int> &indice, Vec<D, int> &numbers,
      int &count, int dim)
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

  template <int D, int ord>
  constexpr Mat<TrefftzElement<D, ord>::npoly, D, int>
  TrefftzElement<D, ord>::MakeIndices ()
  {
    Mat<npoly, D, int> indice = 0;
    Vec<D, int> numbers = 0;
    int count = 0;
    MakeIndices_inner (indice, numbers, count);
    return indice;
  }

  template <int D, int ord>
  constexpr int TrefftzElement<D, ord>::IndexMap (Vec<D, int> index)
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

  template <int D, int ord>
  constexpr Mat<TrefftzElement<D, ord>::nbasis, TrefftzElement<D, ord>::npoly,
                double>
  TrefftzElement<D, ord>::TrefftzBasis ()
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
                i += nbasis - 1;
              }
          }
      }
    return temp_basis;
    // cout << "basis: \n" << basis << endl;
  }

  template <int D, int ord>
  double
  TrefftzElement<D, ord>::ipow_ar (FlatVector<double> base, Vec<D, int> ex,
                                   float result, int count) const
  {
    return count < 0
               ? result
               : ipow_ar (base, ex, pow (base (count), ex (count)) * result,
                          count - 1);
  }

  template <int D, int ord>
  double TrefftzElement<D, ord>::ipowD_ar (int der, FlatVector<double> base,
                                           Vec<D, int> ex, float result,
                                           int count) const
  {
    return ex (der) < 0   ? 0.0
           : count == der ? ipowD_ar (
                 der, base, ex,
                 (ex (count) + 1) * pow (base (count), ex (count)) * result,
                 count - 1)
           : count < 0
               ? result
               : ipowD_ar (der, base, ex,
                           pow (base (count), ex (count)) * result, count - 1);
  }

}

#ifdef NGS_PYTHON
void ExportTrefftzElement (py::module m)
{
  using namespace ngfem;
  py::class_<TrefftzElement<3, 3>, shared_ptr<TrefftzElement<3, 3>>,
             FiniteElement> (m, "TrefftzElement3", "Trefftz space for wave eq")
      .def (py::init<> ())
      .def ("TrefftzBasis", &TrefftzElement<3, 3>::TrefftzBasis)
      //.def("CalcShape", &TrefftzElement<2>::CalcShape)
      .def ("GetNBasis", &TrefftzElement<3, 3>::GetNBasis);
  py::class_<TrefftzElement<2, 3>, shared_ptr<TrefftzElement<2, 3>>,
             FiniteElement> (m, "TrefftzElement2", "Trefftz space for wave eq")
      .def (py::init<> ())
      .def ("TrefftzBasis", &TrefftzElement<2, 3>::TrefftzBasis)
      //.def("CalcShape", &TrefftzElement<2>::CalcShape)
      .def ("GetNBasis", &TrefftzElement<2, 3>::GetNBasis);
}
#endif // NGS_PYTHON
