#include <fem.hpp>
#include "TrefftzElement.hpp"
#include "MultiArray.hpp"
//#include "helpers.cpp"

namespace ngfem
{
  template <int D, int order>
  void TrefftzElement<D, order>::CalcShape (const IntegrationPoint &ip,
                                            BareSliceVector<> shape) const
  {
    for (int l = 0; l < nbasis; l++) // loop over basis functions
      {
        for (int i = 0; i < BinCoeff (D + 1 + order, order);
             i++) // loop over indices
          {
            shape (l) += basisFunctions[l].get (indices[i])
                         * ipow_ar (ip, indices[i]);
          }
      }
  }

  template <int D, int order>
  void TrefftzElement<D, order>::CalcDShape (const IntegrationPoint &ip,
                                             SliceMatrix<> dshape) const
  {
    array<int, D + 1> tempexp;
    int coeff;
    for (int l = 0; l < nbasis; l++) // loop over basis functions
      {
        for (int d = 0; d < D; d++) // loop over derivatives/dimensions
          {
            for (int i = 0; i < BinCoeff (D + 1 + order, order);
                 i++) // loop over indices
              {
                if (indices[i][d + 1] == 0)
                  continue;
                else
                  {
                    tempexp = indices[i];
                    dshape (l, d) += tempexp[d + 1]--
                                     * basisFunctions[l].get (indices[i])
                                     * ipow_ar (ip, tempexp);
                  }
              }
          }
      }
  }

  template <int D, int order>
  void
  TrefftzElement<D, order>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                       BareSliceVector<> shape) const
  {
    for (int l = 0; l < nbasis; l++) // loop over basis functions
      {
        for (int i = 0; i < BinCoeff (D + 1 + order, order);
             i++) // loop over indices
          {
            shape (l) += basisFunctions[l].get (indices[i])
                         * ipow_ar (mip.GetPoint (), indices[i]);
          }
      }
  }

  template <int D, int order>
  void
  TrefftzElement<D, order>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                        SliceMatrix<> dshape) const
  {
    array<int, D + 1> tempexp;
    int coeff;
    for (int l = 0; l < nbasis; l++) // loop over basis functions
      {
        for (int d = 0; d < D; d++) // loop over derivatives/dimensions
          {
            for (int i = 0; i < BinCoeff (D + 1 + order, order);
                 i++) // loop over indices
              {
                if (indices[i][d + 1] == 0)
                  continue;
                else
                  {
                    tempexp = indices[i];
                    dshape (l, d) += tempexp[d + 1]--
                                     * basisFunctions[l].get (indices[i])
                                     * ipow_ar (mip.GetPoint (), tempexp);
                  }
              }
          }
      }
  }

  template <int D, int order> void TrefftzElement<D, order>::TrefftzBasis ()
  {
    for (int l = 0; l < nbasis; l++) // loop over basis functions
      {
        // cout << "======= basis: " << l << endl;

        for (int i = 0; i < BinCoeff (D + 1 + order, order);
             i++) // loop over indices
          {
            float k = (float)indices[i][0];
            if (k > 1)
              {
                // cout << "===== rekursion";
                float temp = 0.0;

                for (int m = 1; m <= D; m++) // rekursive sum
                  {
                    array<int, D + 1> get_coeff = indices[i];
                    get_coeff[0] = get_coeff[0] - 2;
                    get_coeff[m] = get_coeff[m] + 2;
                    temp += (indices[i][m] + 1) * (indices[i][m] + 2)
                            * basisFunctions[l].get (get_coeff);
                  }

                temp = 1 / (k * (k - 1)) * temp;
                basisFunctions[l].put (indices[i], temp);
              }
            else if (k == 0) // time=0
              {
                basisFunctions[l].put (
                    indices[l], 1.0); // set coeff at time=0 to monomial basis
                i += BinCoeff (D + order, order)
                     + BinCoeff (D + order - 1, order - 1);
              }
          }
      }
  }

  template <int D, int order>
  void
  TrefftzElement<D, order>::MakeIndices (int dim, array<int, D + 1> &numbers,
                                         vector<array<int, D + 1>> &indices)
  {
    if (dim > 0)
      {
        for (int i = 0; i <= order; i++)
          {
            numbers[numbers.size () - dim] = i;
            MakeIndices (dim - 1, numbers, indices);
          }
      }
    else
      {
        int sum = 0;
        for (int i = 0; i < numbers.size (); i++)
          {
            sum += numbers[i];
          }
        if (sum <= order)
          {
            indices.push_back (numbers);

            for (int i = 0; i < numbers.size (); i++)
              {
                cout << numbers[i] << " ";
              }
            cout << "\n";
          }
      }
  }

  template <int D, int order>
  float TrefftzElement<D, order>::ipow_ar (IntegrationPoint base,
                                           array<int, D + 1> exp) const
  {
    int result = 1;
    for (int i = 0; i < D + 1; i++)
      {
        result *= pow (base (i), exp[i]);
      }
    return result;
  }

  template <int D, int order>
  float TrefftzElement<D, order>::ipow_ar (FlatVector<double> base,
                                           array<int, D + 1> exp) const
  {
    int result = 1;
    for (int i = 0; i < D + 1; i++)
      {
        result *= pow (base (i), exp[i]);
      }
    return result;
  }

  template <int D, int order> int TrefftzElement<D, order>::GetNBasis ()
  {
    return nbasis;
  }

}

#ifdef NGS_PYTHON
void ExportTrefftzElement (py::module m)
{
  using namespace ngfem;
  py::class_<TrefftzElement<2, 4>, shared_ptr<TrefftzElement<2, 4>>,
             BaseScalarFiniteElement> (m, "TrefftzElement",
                                       "Trefftz space for wave eq")
      //.def(py::init<>())
      .def (py::init<> ())
      .def ("TrefftzBasis", &TrefftzElement<2, 4>::TrefftzBasis)
      //.def("CalcShape", &TrefftzElement<2>::CalcShape)
      .def ("GetNBasis", &TrefftzElement<2, 4>::GetNBasis);
}
#endif // NGS_PYTHON
