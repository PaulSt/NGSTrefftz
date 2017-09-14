#include <fem.hpp>
#include "TrefftzElement.hpp"
//#include "helpers.cpp"

namespace ngfem
{
  template <int D, int ord>
  const Mat<TrefftzElement<D, ord>::npoly, D + 1, int>
      TrefftzElement<D, ord>::indices = MakeIndices ();

  template <int D, int ord>
  const Mat<TrefftzElement<D, ord>::nbasis, TrefftzElement<D, ord>::npoly,
            double>
      TrefftzElement<D, ord>::basis = TrefftzBasis ();

  template <int D, int ord>
  TrefftzElement<D, ord>::TrefftzElement () : FiniteElement ()
  {
    // basisFunctions = vector<MultiArray<float,D+1> >(nbasis, ord+1); //nbasis
    // MultiArrays of depth ord+1 cout << "ord: " + to_string(ord) + ",
    // dimension: " + to_string(D) + ", number of basis functions: " << nbasis
    // << endl; cout << "\n ===== exponentials: \n"; indices.reserve(
    // BinCoeff(D+1+ord,ord) );

    // indices = MakeIndices();
    // basis = TrefftzBasis();
  }

  template <int D, int ord>
  void
  TrefftzElement<D, ord>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                     Vector<double> &shape) const
  {
    Vec<npoly, float> polynomial;

    for (int i = 0; i < npoly; i++) // loop over indices
      {
        polynomial (i) = ipow_ar (mip.GetPoint (), indices.Row (i));
      }

    shape = basis * polynomial;

    // FlatVector<double> point = mip.GetPoint();
    // shape(0) = point(0) * point(1);
  }

  template <int D, int ord>
  void
  TrefftzElement<D, ord>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                      SliceMatrix<> dshape) const
  {
    /*
    array<int, D+1> tempexp;
    int coeff;
    for(int l=0;l<nbasis;l++) //loop over basis functions
    {
            for(int d=0;d<D;d++)  //loop over derivatives/dimensions
            {
                    for(int i=0;i<BinCoeff(D+1 + ord, ord);i++)//loop over
    indices
                    {
                            if(indices[i][d+1] == 0) continue;
                            else
                            {
                                    tempexp = indices[i];
                                    dshape(l,d) += tempexp[d+1]-- *
    basisFunctions[l].get(indices[i]) * ipow_ar(mip.GetPoint(),tempexp,1,D+1);
                            }
                    }
            }
    }
    */
  }

  template <int D, int ord>
  double
  TrefftzElement<D, ord>::ipow_ar (FlatVector<double> base, Vec<D + 1, int> ex,
                                   float result, int count) const
  {
    return count == 0
               ? result
               : ipow_ar (base, ex,
                          pow (base (count - 1), ex (count - 1)) * result,
                          count - 1);
  }

  template <int D, int ord> int TrefftzElement<D, ord>::GetNBasis () const
  {
    return nbasis;
  }

}

#ifdef NGS_PYTHON
void ExportTrefftzElement (py::module m)
{
  using namespace ngfem;
  py::class_<TrefftzElement<2, 3>, shared_ptr<TrefftzElement<2, 3>>,
             FiniteElement> (m, "TrefftzElement", "Trefftz space for wave eq")
      //.def(py::init<>())
      .def (py::init<> ())
      .def ("TrefftzBasis", &TrefftzElement<2, 3>::TrefftzBasis)
      //.def("CalcShape", &TrefftzElement<2>::CalcShape)
      .def ("GetNBasis", &TrefftzElement<2, 3>::GetNBasis);
}
#endif // NGS_PYTHON
