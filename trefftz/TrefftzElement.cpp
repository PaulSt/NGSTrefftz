#include "TrefftzElement.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "recursive_pol.hpp"
#include "helpers.hpp"

#include <ctime>

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
                                         ELEMENT_TYPE aeltype, int basistype)
      : ScalarMappedElement<D> (BinCoeff (D - 1 + aord, aord)
                                    + BinCoeff (D - 1 + aord - 1, aord - 1),
                                aord),
        ord (aord), c (ac), nbasis (BinCoeff (D - 1 + ord, ord)
                                    + BinCoeff (D - 1 + ord - 1, ord - 1)),
        npoly (BinCoeff (D + ord, ord)), indices (GetIndices ()),
        basis (GetTrefftzBasis (basistype)), eltype (aeltype), horner (ord)
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

    // static double elapsed_secs2 = 0;
    // static int calltime2 = 0;
    // clock_t begin2 = clock();
    // 		for(int b = 0; b < nbasis; b++)
    // 			shape(b) = horner.MultiHorner(basis.Row(b),cpoint);
    // clock_t end2 = clock();
    // elapsed_secs2 += double(end2 - begin2) / CLOCKS_PER_SEC;
    // cout << "CS horner: " << elapsed_secs2 <<" times called: " <<
    // ++calltime2 << " avg run: "<< elapsed_secs2/calltime2<<  endl;
  }

  template <int D>
  void T_TrefftzElement<D>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                        SliceMatrix<> dshape) const
  {
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
    dshape.Col (D - 1) *= c;  // inner derivative
    dshape *= (2.0 / elsize); // inner derivative
  }

  template <int D>
  Vector<double> T_TrefftzElement<D>::ShiftPoint (Vector<double> point) const
  {
    point -= elcenter;
    point *= (2.0 / elsize);
    point[D - 1] *= c;
    return point;
  }

  template <int D> Matrix<int> T_TrefftzElement<D>::GetIndices ()
  {
    static int order;
    static Matrix<int> indicesstorage;
    if (order != ord)
      {
        order = ord;
        indicesstorage = MakeIndices ();
      }
    return indicesstorage;
  }

  template <int D>
  Matrix<double> T_TrefftzElement<D>::GetTrefftzBasis (int basistype) const
  {
    static int order;
    static int btype;
    static Matrix<double> basisstorage;
    if (order != ord || btype != basistype)
      {
        basisstorage = TrefftzBasis (basistype);
        order = ord;
        btype = basistype;
      }
    return basisstorage;
  }

  template <int D>
  constexpr Matrix<double>
  T_TrefftzElement<D>::TrefftzBasis (int basistype) const
  {
    Matrix<double> temp_basis (nbasis, npoly);
    temp_basis = 0;
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
                    temp_basis (l, i)
                        += (indices (i, m) + 1) * (indices (i, m) + 2)
                           * temp_basis (l, IndexMap (get_coeff));
                  }
                temp_basis (l, i) *= 1.0 / (k * (k - 1));
              }
            else if (k <= 1 && l == 0
                     && setbasis <= i) // if((k == 0 && l < BinCoeff(D-1 + ord,
                                       // ord)) || (k == 1 && l >= BinCoeff(D-1
                                       // + ord, ord))) //time=0 and =1
              {
                double coeff = 1;
                switch (basistype)
                  {
                  case 0:
                    temp_basis (setbasis++, i) = 1.0; // set the l-th coeff to
                                                      // 1
                    // i += nbasis-1;	//jump to time = 2 if i=0
                    break;
                  case 1:
                    for (int exponent : indices.Row (i).Range (0, D - 1))
                      coeff *= LegCoeffMonBasis (l, exponent);
                    temp_basis (l, i) = coeff;
                    break;
                  case 2:
                    for (int exponent : indices.Row (i).Range (0, D - 1))
                      coeff *= ChebCoeffMonBasis (l, exponent);
                    temp_basis (l, i) = coeff;
                    break;
                  }
              }
          }
      }
    return temp_basis;
  }

  template <int D>
  constexpr void
  T_TrefftzElement<D>::MakeIndices_inner (Matrix<int> &indice,
                                          Vec<D, int> &numbers, int &count,
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

  template <int D> constexpr Matrix<int> T_TrefftzElement<D>::MakeIndices ()
  {
    Matrix<int> indice (npoly, D);
    Vec<D, int> numbers = 0;
    int count = 0;
    for (int o = 0; o <= ord; o++)
      {
        MakeIndices_inner (indice, numbers, count, o);
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
