#include "TrefftzElement.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "recursive_pol.hpp"

namespace ngfem
{
  // template<int D>
  // T_TrefftzElement<D> ::T_TrefftzElement()
  // 	: ScalarMappedElement<D>(D+1, 1),
  // 		ord(1),
  // 		nbasis(D+1),
  // 		npoly(D+1),
  // 		indices(MakeIndices()),
  // 		basis(TrefftzBasis()),
  // 		eltype(ET_TRIG)
  // 		{;}
  // template<int D>
  // T_TrefftzElement<D> ::T_TrefftzElement(int aord, float ac = 1)
  // 	: ScalarMappedElement<D>(BinCoeff(D-1 + aord, aord) + BinCoeff(D-1 +
  // aord-1, aord-1), aord), 		ord(aord), 		c(ac), 		nbasis(BinCoeff(D-1 + ord, ord)
  // + BinCoeff(D-1 + ord-1, ord-1)), 		npoly(BinCoeff(D + ord, ord)),
  // 		indices(MakeIndices()),
  // 		basis(TrefftzBasis()),
  // 		eltype(ET_TRIG)
  // 		{;}

  template <int D>
  T_TrefftzElement<D>::T_TrefftzElement (int aord, float ac,
                                         ELEMENT_TYPE aeltype)
      : ScalarMappedElement<D> (BinCoeff (D - 1 + aord, aord)
                                    + BinCoeff (D - 1 + aord - 1, aord - 1),
                                aord),
        ord (aord), c (ac), nbasis (BinCoeff (D - 1 + ord, ord)
                                    + BinCoeff (D - 1 + ord - 1, ord - 1)),
        npoly (BinCoeff (D + ord, ord)), indices (MakeIndices ()),
        basis (TrefftzBasis ()), eltype (aeltype)
  {
    ;
  }

  template <int D>
  void T_TrefftzElement<D>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                       BareSliceVector<> shape) const
  {
    // Vector<double> cpoint = mip.GetPoint();
    Vector<double> cpoint (D);
    for (int d = 0; d < D; d++)
      cpoint[d] = mip.GetPoint ()[d];
    cpoint = ShiftPoint (cpoint);
    Vector<double> polynomial (npoly);

    for (int i = 0; i < npoly; i++) // loop over indices
      {
        polynomial (i) = ipow_ar (cpoint, indices.Row (i));
      }
    Vector<double> tempshape (nbasis);
    tempshape = basis * polynomial;
    for (int b = 0; b < nbasis; b++)
      shape (b) = tempshape (b);
  }

  template <int D>
  void T_TrefftzElement<D>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                        SliceMatrix<> dshape) const
  {
    // Vector<double> cpoint = mip.GetPoint();
    Vector<double> cpoint (D);
    for (int d = 0; d < D; d++)
      cpoint[d] = mip.GetPoint ()[d];
    cpoint = ShiftPoint (cpoint);
    Vector<double> polynomial (npoly);

    for (int d = 0; d < D; d++) // loop over derivatives/dimensions
      {
        for (int i = 0; i < npoly; i++) // loop over indices
          {
            polynomial (i) = ipowD_ar (d, cpoint, indices.Row (i));
          }
        dshape.Col (d) = basis * polynomial;
      }
    dshape.Col (0) *= c;      // inner derivative
    dshape *= (1.0 / elsize); // inner derivative
  }

  template <int D>
  Vector<double> T_TrefftzElement<D>::ShiftPoint (Vector<double> point) const
  {
    point -= elcenter;
    point *= (1.0 / elsize);
    point[0] *= c;
    return point;
  }

  template <int D>
  constexpr Matrix<double> T_TrefftzElement<D>::TrefftzBasis () const
  {
    Matrix<double> temp_basis (nbasis, npoly);
    temp_basis = 0;
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
            else if (k == 0 && l < BinCoeff (D - 1 + ord, ord)) // time=0 and
                                                                // =1
              {
                // temp_basis( l, IndexMap(indices.Row(l)) ) = 1.0; //set the
                // l-th coeff to 1 i += nbasis-1;	//jump to time = 2
                double coeff = 1;
                for (int exponent : indices.Row (i).Range (1, D))
                  coeff *= LegCoeffMonBasis (l, exponent);
                temp_basis (l, IndexMap (indices.Row (i))) = coeff;
                // LegendrePolynomial leg;
                // cout << "legendre pol: " << endl << leg.GetCoefs() << endl;
              }
            else if (k == 1
                     && l >= BinCoeff (D - 1 + ord, ord)) // time=0 and =1
              {
                double coeff = 1;
                for (int exponent : indices.Row (i).Range (1, D))
                  coeff *= LegCoeffMonBasis (l, exponent);
                temp_basis (l, IndexMap (indices.Row (i))) = coeff;
              }
          }
      }
    return temp_basis;
  }

  template <int D>
  constexpr void T_TrefftzElement<D>::MakeIndices_inner (Matrix<int> &indice,
                                                         Vec<D, int> &numbers,
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
          }
      }
  }

  template <int D> constexpr Matrix<int> T_TrefftzElement<D>::MakeIndices ()
  {
    Matrix<int> indice (npoly, D);
    Vec<D, int> numbers = 0;
    int count = 0;
    MakeIndices_inner (indice, numbers, count);
    return indice;
  }

  template <int D>
  constexpr int T_TrefftzElement<D>::IndexMap (Vec<D, int> index) const
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

  template <int D>
  double T_TrefftzElement<D>::ipow_ar (FlatVector<double> base, Vec<D, int> ex,
                                       double result, int count) const
  {
    return count < 0
               ? result
               : ipow_ar (base, ex, pow (base (count), ex (count)) * result,
                          count - 1);
  }

  template <int D>
  double T_TrefftzElement<D>::ipowD_ar (int der, FlatVector<double> base,
                                        Vec<D, int> ex, double result,
                                        int count) const
  {
    return ex (der) == 0  ? 0.0
           : count == der ? ipowD_ar (
                 der, base, ex,
                 (ex (count)) * pow (base (count), ex (count) - 1) * result,
                 count - 1)
           : count < 0
               ? result
               : ipowD_ar (der, base, ex,
                           pow (base (count), ex (count)) * result, count - 1);
  }

  template class T_TrefftzElement<1>;
  template class T_TrefftzElement<2>;
  template class T_TrefftzElement<3>;
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
