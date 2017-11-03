#include "TrefftzElement.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"

namespace ngfem
{

  void
  BaseScalarMappedElement ::CalcShape (const BaseMappedIntegrationPoint &mip,
                                       BareSliceVector<Complex> shape) const
  {
    CalcShape (mip,
               SliceVector<double> (ndof, 2 * shape.Dist (),
                                    reinterpret_cast<double *> (&shape (0))));
    SliceVector<double> imag_part (
        ndof, 2 * shape.Dist (), reinterpret_cast<double *> (&shape (0)) + 1);
    imag_part = 0.0;
  }

  void
  BaseScalarMappedElement ::CalcShape (const BaseMappedIntegrationRule &mir,
                                       SliceMatrix<> shape) const
  {
    for (int i = 0; i < mir.Size (); i++)
      CalcShape (mir[i], shape.Col (i));
  }

  void BaseScalarMappedElement ::CalcShape (
      const SIMD_BaseMappedIntegrationRule &mir,
      BareSliceMatrix<SIMD<double>> shape) const
  {
    throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
  }

  double
  BaseScalarMappedElement ::Evaluate (const BaseMappedIntegrationPoint &mip,
                                      BareSliceVector<double> x) const
  {
    VectorMem<20, double> shape (ndof);
    CalcShape (mip, shape);
    // cout << "shape: " << shape << endl;
    // cout << "x: " << x(0) << x(1) << endl;
    // cout << "innter prod: " << InnerProduct (shape, x) << endl;
    return InnerProduct (shape, x);
  }

  void
  BaseScalarMappedElement ::Evaluate (const BaseMappedIntegrationRule &mir,
                                      BareSliceVector<double> coefs,
                                      FlatVector<double> vals) const
  {
    for (size_t i = 0; i < mir.Size (); i++) //.GetNIP(); i++)
      vals (i) = Evaluate (mir[i], coefs);
  }

  void BaseScalarMappedElement ::Evaluate (
      const SIMD_BaseMappedIntegrationRule &mir, BareSliceVector<> coefs,
      BareVector<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (
        string ("Evaluate (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void BaseScalarMappedElement ::Evaluate (
      const SIMD_BaseMappedIntegrationRule &mir, SliceMatrix<> coefs,
      BareSliceMatrix<SIMD<double>> values) const
  {
    for (size_t i = 0; i < coefs.Width (); i++)
      Evaluate (mir, coefs.Col (i), values.Row (i));
  }

  void
  BaseScalarMappedElement ::Evaluate (const BaseMappedIntegrationRule &mir,
                                      SliceMatrix<> coefs,
                                      SliceMatrix<> values) const
  {
    VectorMem<100> shapes (coefs.Height ());
    for (size_t i = 0; i < mir.Size (); i++)
      {
        CalcShape (mir[i], shapes);
        values.Row (i) = Trans (coefs) * shapes;
      }
  }

  void BaseScalarMappedElement ::EvaluateTrans (
      const BaseMappedIntegrationRule &mir, FlatVector<double> vals,
      BareSliceVector<double> coefs) const
  {
    VectorMem<20, double> shape (ndof);
    coefs.AddSize (ndof) = 0.0;
    for (int i = 0; i < mir.Size (); i++) // GetNIP()
      {
        CalcShape (mir[i], shape);
        coefs.AddSize (ndof) += vals (i) * shape;
      }
  }

  void BaseScalarMappedElement ::AddTrans (
      const SIMD_BaseMappedIntegrationRule &mir,
      BareVector<SIMD<double>> values, BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD (
        string ("AddTrans (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void BaseScalarMappedElement ::AddTrans (
      const SIMD_BaseMappedIntegrationRule &mir,
      BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const
  {
    for (int i = 0; i < coefs.Width (); i++)
      AddTrans (mir, values.Row (i), coefs.Col (i));
  }

  void BaseScalarMappedElement ::CalcMappedDShape (
      const SIMD_BaseMappedIntegrationRule &mir,
      BareSliceMatrix<SIMD<double>> dshapes) const
  {
    throw ExceptionNOSIMD ("SIMD - CalcDShape not overloaded");
  }

  void BaseScalarMappedElement ::EvaluateGrad (
      const SIMD_BaseMappedIntegrationRule &ir, BareSliceVector<> coefs,
      BareSliceMatrix<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (
        string ("EvaluateGrad (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void BaseScalarMappedElement ::EvaluateGrad (
      const SIMD_IntegrationRule &ir, BareSliceVector<> coefs,
      BareSliceMatrix<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (
        string ("EvaluateGrad (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void BaseScalarMappedElement ::AddGradTrans (
      const SIMD_BaseMappedIntegrationRule &ir,
      BareSliceMatrix<SIMD<double>> values, BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD (
        string ("AddGradTrans (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D> string ScalarMappedElement<D>::ClassName () const
  {
    return "ScalarMappedElement";
  }

  template <int D>
  void ScalarMappedElement<D>::CalcMappedDShape (
      const MappedIntegrationPoint<D, D> &mip, SliceMatrix<> dshape) const
  {
    CalcDShape (mip, dshape);
    /* no mapping - no inner derivative
for (int i = 0; i < dshape.Height(); i++)
{
Vec<D> hv = dshape.Row(i);
FlatVec<D> (&dshape(i,0)) = Trans (mip.GetJacobianInverse ()) * hv;
}
    */
  }

  template <int D>
  void ScalarMappedElement<D>::CalcMappedDShape (
      const MappedIntegrationRule<D, D> &mir, SliceMatrix<> dshapes) const
  {
    for (int i = 0; i < mir.Size (); i++)
      CalcMappedDShape (mir[i], dshapes.Cols (i * D, (i + 1) * D));
  }

  template <int D>
  Vec<D>
  ScalarMappedElement<D>::EvaluateGrad (const BaseMappedIntegrationPoint &ip,
                                        BareSliceVector<double> x) const
  {
    MatrixFixWidth<D> dshape (ndof);
    CalcDShape (ip, dshape);
    Vec<D> grad = Trans (dshape) * x;
    return grad;
  }

  template <int D>
  void ScalarMappedElement<D>::EvaluateGrad (
      const BaseMappedIntegrationRule &ir, BareSliceVector<double> coefs,
      FlatMatrixFixWidth<D, double> vals) const
  {
    for (size_t i = 0; i < ir.Size (); i++)
      vals.Row (i) = EvaluateGrad (ir[i], coefs);
  }

  template <int D>
  void ScalarMappedElement<D>::EvaluateGradTrans (
      const BaseMappedIntegrationRule &ir, FlatMatrixFixWidth<D, double> vals,
      BareSliceVector<double> coefs) const
  {
    MatrixFixWidth<D> dshape (ndof);
    coefs.AddSize (ndof) = 0.0;
    for (int i = 0; i < ir.Size (); i++)
      {
        CalcDShape (ir[i], dshape);
        coefs.AddSize (ndof) += dshape * vals.Row (i);
      }
  }

  template <int D>
  void ScalarMappedElement<D>::EvaluateGradTrans (
      const BaseMappedIntegrationRule &ir, SliceMatrix<> values,
      SliceMatrix<> coefs) const
  {
#ifndef __CUDA_ARCH__
    cout << "EvalGradTrans not overloaded" << endl;
#endif
  }

  template <int D>
  void
  ScalarMappedElement<D>::GetPolOrders (FlatArray<PolOrder<D>> orders) const
  {
#ifndef __CUDA_ARCH__
    throw Exception (string ("GetPolOrders not implemnted for element")
                     + ClassName ());
#endif
  }

  /*
                  template<int D>
                  void ScalarMappedElement<D> :: CalcDDShape (const
     IntegrationPoint & ip, FlatMatrix<> ddshape) const
                  {
                          int nd = GetNDof();
                          int sdim = D;

                          double eps = 1e-7;
                          Matrix<> dshape1(nd, sdim), dshape2(nd, sdim);

                          for (int i = 0; i < sdim; i++)
                                  {
                  IntegrationPoint ip1 = ip;
                  IntegrationPoint ip2 = ip;

                                          ip1(i) -= eps;
                                          ip2(i) += eps;

                  CalcDShape (ip1, dshape1);
                  CalcDShape (ip2, dshape2);
                  dshape2 -= dshape1;
                  dshape2 *= (0.5 / eps);
                  for (int j = 0; j < nd; j++)
                          for (int k = 0; k < sdim; k++)
                                  ddshape(j,sdim*i+k) = dshape2(j,k);
                                  }
                  }


                  template<int D>
                  void ScalarMappedElement<D> :: CalcMappedDDShape (const
     MappedIntegrationPoint<D,D> & mip, SliceMatrix<> ddshape) const
                  {
                          int nd = GetNDof();

                          double eps = 1e-7;
                          MatrixFixWidth<D> dshape1(nd), dshape2(nd);
                          const ElementTransformation & eltrans =
     mip.GetTransformation();

                          for (int i = 0; i < D; i++)
                                  {
                  IntegrationPoint ip1 = mip.IP();
                  IntegrationPoint ip2 = mip.IP();
                                          ip1(i) -= eps;
                                          ip2(i) += eps;
                                          MappedIntegrationPoint<D,D> mip1(ip1,
     eltrans); MappedIntegrationPoint<D,D> mip2(ip2, eltrans);

                  CalcMappedDShape (mip1, dshape1);
                  CalcMappedDShape (mip2, dshape2);

                                          ddshape.Cols(D*i,D*(i+1)) = (0.5/eps)
     * (dshape2-dshape1);
                                  }

                          for (int j = 0; j < D; j++)
                                  {
                                          for (int k = 0; k < nd; k++)
                                                  for (int l = 0; l < D; l++)
                                                          dshape1(k,l) =
     ddshape(k, l*D+j);

                                          dshape2 = dshape1 *
     mip.GetJacobianInverse();

                                          for (int k = 0; k < nd; k++)
                                                  for (int l = 0; l < D; l++)
                                                          ddshape(k, l*D+j) =
     dshape2(k,l);
                                  }

                  }
          **/

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    FlatVector<double> point = ShiftPoint (mip.GetPoint ());
    Vec<npoly, float> polynomial;

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
    FlatVector<double> point = ShiftPoint (mip.GetPoint ());
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
        if (d == 0)
          dshape.Col (d) *= c; // inner derivative
      }
    dshape *= (1.0 / elsize); // inner derivative
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D, int ord>
  const Matrix<double>
      TrefftzHelmholtzElement<D, ord>::directions = MakeDirections ();

  template <int D, int ord>
  void TrefftzHelmholtzElement<D, ord>::CalcShape (
      const BaseMappedIntegrationPoint &mip, BareSliceVector<> shape) const
  {
    dcomp i = sqrt (-1);
    for (int i = 0; i < nbasis; i++)
      {
        cout << exp (i * InnerProduct (directions.Row (i), mip.GetPoint ()));
      }
  }

  template <int D, int ord>
  void TrefftzHelmholtzElement<D, ord>::CalcDShape (
      const BaseMappedIntegrationPoint &mip, SliceMatrix<> dshape) const
  {
  }

  template <int D, int ord>
  constexpr Mat<TrefftzHelmholtzElement<D, ord>::nbasis, D, double>
  TrefftzHelmholtzElement<D, ord>::MakeDirections ()
  {
    Mat<nbasis, D, double> dir = 0;
    float theta = 0;
    for (int i = 0; i < nbasis; i++)
      {
        theta = (2 * M_PI / nbasis) * i;
        dir (i, 0) = cos (theta);
        dir (i, 1) = sin (theta);
      }
    return dir;
    // cout << "basis: \n" << basis << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NGS_PYTHON
void ExportTrefftzElement (py::module m)
{
  using namespace ngfem;
  py::class_<TrefftzElement<3, 3>, shared_ptr<TrefftzElement<3, 3>>,
             FiniteElement> (m, "TrefftzElement3", "Trefftz space for wave eq")
      .def (py::init<> ());
  py::class_<TrefftzElement<2, 3>, shared_ptr<TrefftzElement<2, 3>>,
             FiniteElement> (m, "TrefftzElement2", "Trefftz space for wave eq")
      .def (py::init<> ());

  py::class_<TrefftzHelmholtzElement<2, 3>,
             shared_ptr<TrefftzHelmholtzElement<2, 3>>, FiniteElement> (
      m, "TrefftzHelmholtzElement2", "Trefftz space for Helmholtz eq")
      .def (py::init<> ());
}
#endif // NGS_PYTHON
