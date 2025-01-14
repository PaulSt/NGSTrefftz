#include "scalarmappedfe.hpp"

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
                                       BareSliceMatrix<> shape) const
  {
    for (size_t i = 0; i < mir.Size (); i++)
      CalcShape (mir[i], shape.Col (i));
  }

  void
  BaseScalarMappedElement ::CalcShape (const SIMD_BaseMappedIntegrationRule &,
                                       BareSliceMatrix<SIMD<double>>) const
  {
    cout << "SIMD - CalcShape not overloaded" << endl;
    throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
  }

  double
  BaseScalarMappedElement ::Evaluate (const BaseMappedIntegrationPoint &mip,
                                      BareSliceVector<double> x) const
  {
    VectorMem<20, double> shape (ndof);
    CalcShape (mip, shape);
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

  void
  BaseScalarMappedElement ::Evaluate (const SIMD_BaseMappedIntegrationRule &,
                                      BareSliceVector<>,
                                      BareVector<SIMD<double>>) const
  {
    cout << "SIMD - Eval not overloaded" << endl;
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
    coefs.Range (0, ndof) = 0.0;
    for (size_t i = 0; i < mir.Size (); i++) // GetNIP()
      {
        CalcShape (mir[i], shape);
        coefs.Range (0, ndof) += vals (i) * shape;
      }
  }

  void
  BaseScalarMappedElement ::AddTrans (const SIMD_BaseMappedIntegrationRule &,
                                      BareVector<SIMD<double>>,
                                      BareSliceVector<>) const
  {
    cout << "SIMD - AddTrans not overloaded" << endl;
    throw ExceptionNOSIMD (
        string ("AddTrans (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void BaseScalarMappedElement ::AddTrans (
      const SIMD_BaseMappedIntegrationRule &mir,
      BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const
  {
    for (size_t i = 0; i < coefs.Width (); i++)
      AddTrans (mir, values.Row (i), coefs.Col (i));
  }

  void BaseScalarMappedElement ::CalcMappedDShape (
      const SIMD_BaseMappedIntegrationRule &,
      BareSliceMatrix<SIMD<double>>) const
  {
    cout << "SIMD - CalcMappedDShape not overloaded" << endl;
    throw ExceptionNOSIMD ("SIMD - CalcDShape not overloaded");
  }

  void BaseScalarMappedElement ::EvaluateGrad (
      const SIMD_BaseMappedIntegrationRule &, BareSliceVector<>,
      BareSliceMatrix<SIMD<double>>) const
  {
    cout << "SIMD - EvaluateGrad not overloaded" << endl;
    throw ExceptionNOSIMD (
        string ("EvaluateGrad (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  // void BaseScalarMappedElement ::
  // EvaluateGrad (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs,
  // BareSliceMatrix<SIMD<double>> values) const
  //{
  // throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for
  // class ")+typeid(*this).name());
  // }

  void BaseScalarMappedElement ::AddGradTrans (
      const SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<double>>,
      BareSliceVector<>) const
  {
    cout << "SIMD - AddTransGrad not overloaded" << endl;
    throw ExceptionNOSIMD (
        string ("AddGradTrans (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void MatToCSR (Matrix<> mat, CSR &sparsemat)
  {
    int spsize = 0;
    sparsemat[0].SetSize (0);
    sparsemat[1].SetSize (0);
    sparsemat[2].SetSize (0);
    for (size_t i = 0; i < mat.Height (); i++)
      {
        // guarantee sparsemat[0].Size() to be nbasis + 1,
        // could be shorter but easier for mat*vec
        sparsemat[0].Append (spsize);
        for (size_t j = 0; j < mat.Width (); j++)
          {
            if (mat (i, j))
              {
                spsize++;
                sparsemat[1].Append (j);
                sparsemat[2].Append (mat (i, j));
              }
          }
      }
    sparsemat[0].Append (spsize);
  };

  template <int D> string ScalarMappedElement<D>::ClassName () const
  {
    return "ScalarMappedElement";
  }

  template <int D>
  void
  ScalarMappedElement<D>::CalcDShape (const BaseMappedIntegrationRule &mir,
                                      BareSliceMatrix<> dshapes) const
  {
    for (size_t i = 0; i < mir.Size (); i++)
      CalcDShape (mir[i], dshapes.Cols (i * D, (i + 1) * D));
  }

  // template<int D>
  // void ScalarMappedElement<D> ::
  // CalcMappedDShape (const BaseMappedIntegrationPoint & mip,
  // BareSliceMatrix<> dshape) const
  //{
  ////no mapping - no inner derivative
  // CalcDShape (mip, dshape);
  //}

  template <int D>
  void ScalarMappedElement<D>::CalcMappedDShape (
      const BaseMappedIntegrationRule &mir, BareSliceMatrix<> dshapes) const
  {
    for (size_t i = 0; i < mir.Size (); i++)
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
    coefs.Range (0, ndof) = 0.0;
    for (size_t i = 0; i < ir.Size (); i++)
      {
        CalcDShape (ir[i], dshape);
        coefs.Range (0, ndof) += dshape * vals.Row (i);
      }
  }

  template <int D>
  void
  ScalarMappedElement<D>::EvaluateGradTrans (const BaseMappedIntegrationRule &,
                                             SliceMatrix<>,
                                             SliceMatrix<>) const
  {
#ifndef __CUDA_ARCH__
    cout << "EvalGradTrans not overloaded" << endl;
#endif
  }

  // template<int D>
  // void ScalarMappedElement<D> ::
  // GetPolOrders (FlatArray<PolOrder<D> > orders) const
  //{
  // #ifndef __CUDA_ARCH__
  // throw Exception (string ("GetPolOrders not implemnted for element") +
  // ClassName());
  // #endif
  //}

  // template<int D>
  // void ScalarMappedElement<D> ::
  // CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
  // BareSliceMatrix<SIMD<double>> dshape) const
  //{
  // cout<<"SIMD - CalcDShape not overloaded"<< endl;
  // throw ExceptionNOSIMD("SIMD - CalcDShape not overloaded");
  //}

  // template<int D>
  // void ScalarMappedElement<D> :: CalcMappedDDShape (const
  // BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> hddshape) const
  //{
  // auto & mip = static_cast<const MappedIntegrationPoint<D,D> &> (bmip);
  // int nd = GetNDof();
  // auto ddshape = hddshape.AddSize(nd, D*D);
  // double eps = 1e-7;
  // MatrixFixWidth<(D)> dshape1(nd), dshape2(nd);
  // const ElementTransformation & eltrans = mip.GetTransformation();

  // for (int i = 0; i < D; i++)
  //{
  // IntegrationPoint ip1 = mip.IP();
  // IntegrationPoint ip2 = mip.IP();
  // ip1(i) -= eps;
  // ip2(i) += eps;
  // MappedIntegrationPoint<D,D> mip1(ip1, eltrans);
  // MappedIntegrationPoint<D,D> mip2(ip2, eltrans);

  ////cout << bmip.GetPoint() << endl;
  ////bmip.GetPoint()(i) -= eps;
  ////cout << bmip.GetPoint() << endl;
  // CalcMappedDShape (mip1, dshape1);
  ////bmip.GetPoint()(i) += 2*eps;
  // CalcMappedDShape (mip2, dshape2);
  ////bmip.GetPoint()(i) -= eps;

  // ddshape.Cols(D*i,D*(i+1)) = (0.5/eps) * (dshape2-dshape1);
  //}

  // for (int j = 0; j < D; j++)
  //{
  // for (int k = 0; k < nd; k++)
  // for (int l = 0; l < D; l++)
  // dshape1(k,l) = ddshape(k, l*D+j);

  // dshape2 = dshape1 * mip.GetJacobianInverse();

  // for (int k = 0; k < nd; k++)
  // for (int l = 0; l < D; l++)
  // ddshape(k, l*D+j) = dshape2(k,l);
  //}
  //}

  template <>
  void ScalarMappedElement<1>::CalcMappedDDShape (
      const BaseMappedIntegrationPoint &, BareSliceMatrix<>) const
  {
    ;
  }

  template <>
  void ScalarMappedElement<2>::CalcMappedDDShape (
      const BaseMappedIntegrationPoint &bmip, BareSliceMatrix<> hddshape) const
  {
    auto ddshape = hddshape.AddSize (this->ndof, 2 * 2);

    // auto & mip = static_cast<const MappedIntegrationPoint<2,2> &> (bmip);
    Vec<2> cpoint = bmip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 2 * (order + 1) + 1);
    mem2[0] = 0;
    double *polxt[2];
    for (int d = 0; d < 2; d++)
      {
        polxt[d] = &mem2[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d1 = 0; d1 < 2; d1++)
      {
        for (int d2 = 0; d2 < 2; d2++)
          {
            Vector<double> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                pol[ii++] = (d1 != d2 ? i * j
                                      : (d1 == 0 ? i * (i - 1) : j * (j - 1)))
                            * polxt[0][i - (d1 == 0) - (d2 == 0)]
                            * polxt[1][j - (d1 == 1) - (d2 == 1)];

            for (int i = 0; i < this->ndof; ++i)
              {
                ddshape (i, d2 * 2 + d1) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  ddshape (i, d2 * 2 + d1) += (localmat)[2][j]
                                              * pol[(localmat)[1][j]]
                                              * scale[d1] * scale[d2];
              }
          }
      }
  }

  template <>
  void ScalarMappedElement<3>::CalcMappedDDShape (
      const BaseMappedIntegrationPoint &bmip, BareSliceMatrix<> hddshape) const
  {
    auto ddshape = hddshape.AddSize (this->ndof, 3 * 3);
    // auto & mip = static_cast<const MappedIntegrationPoint<2,2> &> (bmip);
    Vec<3> cpoint = bmip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 3 * (order + 1) + 1);
    mem2[0] = 0;
    double *polxt[3];
    for (int d = 0; d < 3; d++)
      {
        polxt[d] = &mem2[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d1 = 0; d1 < 3; d1++)
      {
        for (int d2 = 0; d2 < 3; d2++)
          {
            Vector<double> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                for (int k = 0; k <= order - i - j; k++)
                  {
                    int ijkd1 = d1 == 0 ? i : (d1 == 1 ? j : k);
                    int ijkd2 = d2 == 0 ? i : (d2 == 1 ? j : k);
                    pol[ii++] = ijkd1 * (ijkd2 - (d1 == d2))
                                * polxt[0][i - (d1 == 0) - (d2 == 0)]
                                * polxt[1][j - (d1 == 1) - (d2 == 1)]
                                * polxt[2][k - (d1 == 2) - (d2 == 2)];
                  }

            for (int i = 0; i < this->ndof; ++i)
              {
                ddshape (i, d2 * 3 + d1) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  ddshape (i, d2 * 3 + d1) += (localmat)[2][j]
                                              * pol[(localmat)[1][j]]
                                              * scale[d1] * scale[d2];
              }
          }
      }
  }

  template <>
  void ScalarMappedElement<4>::CalcMappedDDShape (
      const BaseMappedIntegrationPoint &, BareSliceMatrix<>) const
  {
    ;
  }

  ////////////////////////////////////////////////////////
  /////////////// CalcShape implementation ///////////////
  ////////////////////////////////////////////////////////

  template <>
  void
  ScalarMappedElement<1>::CalcShape (const SIMD_BaseMappedIntegrationRule &,
                                     BareSliceMatrix<SIMD<double>>) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void ScalarMappedElement<2>::CalcShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> shape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);
        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 2 * (order + 1));
        Vec<2, SIMD<double> *> polxt;
        for (int d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (order + 1)];
            Monomial (order, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            pol[ii++] = polxt[0][i] * polxt[1][j];
        // TB*monomials for trefftz shape fcts
        for (int i = 0; i < this->ndof; ++i)
          {
            shape (i, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (i, imip) += (localmat)[2][j] * pol[(localmat)[1][j]];
          }
      }
  }

  template <>
  void ScalarMappedElement<3>::CalcShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> shape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);
        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 3 * (order + 1));
        Vec<3, SIMD<double> *> polxt;
        for (int d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (order + 1)];
            Monomial (order, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];
        // TB*monomials for trefftz shape fcts
        for (int i = 0; i < this->ndof; ++i)
          {
            shape (i, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (i, imip) += (localmat)[2][j] * pol[(localmat)[1][j]];
          }
      }
  }

  template <>
  void ScalarMappedElement<4>::CalcShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> shape) const
  {
    // auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&>
    // (mir);
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<4, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);
        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 4 * (order + 1));
        Vec<4, SIMD<double> *> polxt;
        for (int d = 0; d < 4; d++)
          {
            polxt[d] = &mem[d * (order + 1)];
            Monomial (order, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              for (int l = 0; l <= order - i - j - k; l++)
                pol[ii++]
                    = polxt[0][i] * polxt[1][j] * polxt[2][k] * polxt[3][l];
        // TB*monomials for trefftz shape fcts
        for (int i = 0; i < this->ndof; ++i)
          {
            shape (i, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (i, imip) += (localmat)[2][j] * pol[(localmat)[1][j]];
          }
      }
  }

  template <>
  void
  ScalarMappedElement<1>::CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                                      BareSliceMatrix<SIMD<double>>) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void ScalarMappedElement<2>::CalcDShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);
        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 2 * (order + 1) + 1);
        mem[0] = 0;
        Vec<2, SIMD<double> *> polxt;
        for (int d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 2; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                            * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

            for (int i = 0; i < this->ndof; ++i)
              {
                dshape (i * 2 + d, imip) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (i * 2 + d, imip)
                      += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
              }
          }
      }
  }

  template <>
  void ScalarMappedElement<3>::CalcDShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 3 * (order + 1) + 1);
        mem[0] = 0;
        Vec<3, SIMD<double> *> polxt;
        for (int d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 3; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                for (int k = 0; k <= order - i - j; k++)
                  pol[ii++] = (d == 0 ? i : (d == 1 ? j : (d == 2 ? k : 0)))
                              * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                              * polxt[2][k - (d == 2)];

            for (int i = 0; i < this->ndof; ++i)
              {
                dshape (i * 3 + d, imip) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (i * 3 + d, imip)
                      += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
              }
          }
      }
  }

  template <>
  void ScalarMappedElement<4>::CalcDShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<4, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);
        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 4 * (order + 1) + 1);
        mem[0] = 0;
        Vec<4, SIMD<double> *> polxt;
        for (int d = 0; d < 4; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 4; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                for (int k = 0; k <= order - i - j; k++)
                  for (int l = 0; l <= order - i - j - k; l++)
                    pol[ii++]
                        = (d == 0 ? i
                                  : (d == 1 ? j
                                            : (d == 2 ? k : (d == 3 ? l : 0))))
                          * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                          * polxt[2][k - (d == 2)] * polxt[3][l - (d == 3)];

            for (int i = 0; i < this->ndof; ++i)
              {
                dshape (i * 4 + d, imip) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (i * 4 + d, imip)
                      += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
              }
          }
      }
  }

  /////////////// non-simd

  template <>
  void ScalarMappedElement<1>::CalcShape (const BaseMappedIntegrationPoint &,
                                          BareSliceVector<>) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void
  ScalarMappedElement<2>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceVector<> shape) const
  {
    // auto & smir = static_cast<const SIMD_BaseMappedIntegrationRuleD>&>
    // (mir);
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 2 * (order + 1));
    double *polxt[2];
    for (int d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (order + 1)];
        Monomial (order, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (int i = 0, ii = 0; i <= order; i++)
      for (int j = 0; j <= order - i; j++)
        pol[ii++] = polxt[0][i] * polxt[1][j];
    // TB*monomials for trefftz shape fcts
    for (int i = 0; i < this->ndof; ++i)
      {
        shape (i) = 0.0;
        for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
          shape (i) += (localmat)[2][j] * pol[(localmat)[1][j]];
      }
  }

  template <>
  void
  ScalarMappedElement<3>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceVector<> shape) const
  {
    Vec<3> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 3 * (order + 1));
    // double* polxt[3];
    Vec<3, double *> polxt;
    for (int d = 0; d < 3; d++)
      {
        polxt[d] = &mem[d * (order + 1)];
        Monomial (order, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (int i = 0, ii = 0; i <= order; i++)
      for (int j = 0; j <= order - i; j++)
        for (int k = 0; k <= order - i - j; k++)
          pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];
    // TB*monomials for trefftz shape fcts
    for (int i = 0; i < this->ndof; ++i)
      {
        shape (i) = 0.0;
        for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
          shape (i) += (localmat)[2][j] * pol[(localmat)[1][j]];
      }
  }

  template <>
  void
  ScalarMappedElement<4>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceVector<> shape) const
  {
    Vec<4> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 4 * (order + 1));
    double *polxt[4];
    for (int d = 0; d < 4; d++)
      {
        polxt[d] = &mem[d * (order + 1)];
        Monomial (order, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (int i = 0, ii = 0; i <= order; i++)
      for (int j = 0; j <= order - i; j++)
        for (int k = 0; k <= order - i - j; k++)
          for (int l = 0; l <= order - i - j - k; l++)
            pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k] * polxt[3][l];
    // TB*monomials for trefftz shape fcts
    for (int i = 0; i < this->ndof; ++i)
      {
        shape (i) = 0.0;
        for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
          shape (i) += (localmat)[2][j] * pol[(localmat)[1][j]];
      }
  }

  template <>
  void ScalarMappedElement<1>::CalcDShape (const BaseMappedIntegrationPoint &,
                                           BareSliceMatrix<>) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void
  ScalarMappedElement<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                      BareSliceMatrix<> dshape) const
  {
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 2 * (order + 1) + 1);
    mem[0] = 0;
    double *polxt[2];
    for (int d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 2; d++)
      {
        Vector<double> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                        * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i, d)
                  += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
          }
      }
  }

  template <>
  void
  ScalarMappedElement<3>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                      BareSliceMatrix<> dshape) const
  {
    Vec<3> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 3 * (order + 1) + 1);
    mem[0] = 0;
    // double* polxt[4];
    Vec<3, double *> polxt;
    for (int d = 0; d < 3; d++)
      {
        polxt[d] = &mem[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 3; d++)
      {
        Vector<double> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              pol[ii++] = (d == 0 ? i : (d == 1 ? j : (d == 2 ? k : 0)))
                          * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                          * polxt[2][k - (d == 2)];

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i, d)
                  += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
          }
      }
  }

  template <>
  void
  ScalarMappedElement<4>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                      BareSliceMatrix<> dshape) const
  {
    Vec<4> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 4 * (order + 1) + 1);
    mem[0] = 0;
    double *polxt[4];
    for (int d = 0; d < 4; d++)
      {
        polxt[d] = &mem[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 4; d++)
      {
        Vector<double> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              for (int l = 0; l <= order - i - j - k; l++)
                pol[ii++]
                    = (d == 0 ? i
                              : (d == 1 ? j : (d == 2 ? k : (d == 3 ? l : 0))))
                      * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                      * polxt[2][k - (d == 2)] * polxt[3][l - (d == 3)];

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i, d)
                  += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
          }
      }
  }

  template <>
  void ScalarMappedElement<1>::CalcDDWaveOperator (
      const SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<double>>,
      BareSliceMatrix<SIMD<double>>, BareSliceMatrix<SIMD<double>>) const
  {
    ;
  }
  template <>
  void ScalarMappedElement<4>::CalcDDWaveOperator (
      const SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<double>>,
      BareSliceMatrix<SIMD<double>>, BareSliceMatrix<SIMD<double>>) const
  {
    ;
  }

  template <>
  void ScalarMappedElement<2>::CalcDDWaveOperator (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape,
      BareSliceMatrix<SIMD<double>> wavespeed,
      BareSliceMatrix<SIMD<double>> mu) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);

        STACK_ARRAY (SIMD<double>, mem, 2 * (order + 1) + 2);
        mem[0] = 0;
        mem[1] = 0;
        Vec<2, SIMD<double> *> polxt;
        for (int d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 2];
            Monomial (order, cpoint[d], polxt[d]);
          }

        Vector<SIMD<double>> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            pol[ii++] = (i * (i - 1) * polxt[0][i - 2] * polxt[1][j]
                             * pow (scale[0], 2)
                         - j * (j - 1) * polxt[0][i] * polxt[1][j - 2]
                               * pow (scale[1], 2) * wavespeed (0, imip))
                        * mu (0, imip);

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i * 2, imip) = 0.0;
            dshape (i * 2 + 1, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i * 2 + 1, imip)
                  += (localmat)[2][j] * pol[(localmat)[1][j]];
          }
      }
  }

  template <>
  void ScalarMappedElement<3>::CalcDDWaveOperator (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape,
      BareSliceMatrix<SIMD<double>> wavespeed,
      BareSliceMatrix<SIMD<double>> mu) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);

        STACK_ARRAY (SIMD<double>, mem, 3 * (order + 1) + 2);
        mem[0] = 0;
        mem[1] = 0;
        Vec<3, SIMD<double> *> polxt;
        for (int d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 2];
            Monomial (order, cpoint[d], polxt[d]);
          }

        Vector<SIMD<double>> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              pol[ii++] = (i * (i - 1) * polxt[0][i - 2] * polxt[1][j]
                               * polxt[2][k] * pow (scale[0], 2)
                           + j * (j - 1) * polxt[0][i] * polxt[1][j - 2]
                                 * polxt[2][k] * pow (scale[1], 2)
                           - k * (k - 1) * polxt[0][i] * polxt[1][j]
                                 * polxt[2][k - 2] * pow (scale[2], 2)
                                 * wavespeed (0, imip))
                          * mu (0, imip);

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i * 3, imip) = 0.0;
            dshape (i * 3 + 1, imip) = 0.0;
            dshape (i * 3 + 2, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i * 3 + 2, imip)
                  += (localmat)[2][j] * pol[(localmat)[1][j]];
          }
      }
  }

  template class ScalarMappedElement<1>;
  template class ScalarMappedElement<2>;
  template class ScalarMappedElement<3>;
  template class ScalarMappedElement<4>;

  template <int D>
  void
  BlockMappedElement<D>::CalcShape (const SIMD_BaseMappedIntegrationRule &,
                                    BareSliceMatrix<SIMD<double>>) const
  {
    throw Exception ("Not allowed");
  }

  template <int D>
  void BlockMappedElement<D>::CalcShape (const BaseMappedIntegrationPoint &,
                                         BareSliceVector<>) const
  {
    throw Exception ("Not allowed");
  }

  template <>
  void BlockMappedElement<1>::CalcDShape (const BaseMappedIntegrationPoint &,
                                          BareSliceMatrix<>) const
  {
    throw Exception ("Not implemented point ");
  }
  template <>
  void BlockMappedElement<4>::CalcDShape (const BaseMappedIntegrationPoint &,
                                          BareSliceMatrix<>) const
  {
    throw Exception ("Not implemented point ");
  }

  template <>
  void
  BlockMappedElement<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceMatrix<> dshape) const
  {
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 2 * (order + 1) + 1);
    mem[0] = 0;
    double *polxt[2];
    for (int d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 2; d++)
      {
        Vector<double> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            pol[ii++] = polxt[0][i] * polxt[1][j];

        CSR localmat = localmats[d];
        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i, d)
                  += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
          }
      }
  }

  template <>
  void
  BlockMappedElement<3>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceMatrix<> dshape) const
  {
    Vec<3> cpoint = mip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 3 * (order + 1) + 1);
    mem[0] = 0;
    // double* polxt[4];
    Vec<3, double *> polxt;
    for (int d = 0; d < 3; d++)
      {
        polxt[d] = &mem[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 3; d++)
      {
        Vector<double> pol (npoly);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];

        CSR localmat = localmats[d];
        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i, d)
                  += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
          }
      }
  }

  // template<int D>
  // void BlockMappedElement<D> :: CalcDShape (const BaseMappedIntegrationRule
  // & mir, BareSliceMatrix<> dshapes) const
  //{ throw Exception("Not implemented no simd "); }

  template <>
  void
  BlockMappedElement<1>::CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                                     BareSliceMatrix<SIMD<double>>) const
  {
    throw Exception ("Not implemented 1d");
  }

  template <>
  void BlockMappedElement<2>::CalcDShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);
        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 2 * (order + 1) + 1);
        mem[0] = 0;
        Vec<2, SIMD<double> *> polxt;
        for (int d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 2; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                pol[ii++] = polxt[0][i] * polxt[1][j];

            CSR localmat = localmats[d];
            for (int i = 0; i < ndof; ++i)
              {
                dshape (i * 2 + d, imip) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (i * 2 + d, imip)
                      += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
              }
          }
      }
  }

  template <>
  void BlockMappedElement<3>::CalcDShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 3 * (order + 1) + 1);
        mem[0] = 0;
        Vec<3, SIMD<double> *> polxt;
        for (int d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 3; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                for (int k = 0; k <= order - i - j; k++)
                  pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];

            CSR localmat = localmats[d];
            for (int i = 0; i < ndof; ++i)
              {
                dshape (i * 3 + d, imip) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (i * 3 + d, imip)
                      += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
              }
          }
      }
  }

  template <>
  void BlockMappedElement<4>::CalcDShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (size_t imip = 0; imip < smir.Size (); imip++)
      {
        Vec<4, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= shift;
        vtimes (cpoint, scale);

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 4 * (order + 1) + 1);
        mem[0] = 0;
        Vec<4, SIMD<double> *> polxt;
        for (size_t d = 0; d < 4; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 4; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                for (int k = 0; k <= order - i - j; k++)
                  for (int l = 0; l <= order - i - j - k; l++)
                    pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k]
                                * polxt[3][l];

            CSR localmat = localmats[d];
            for (int i = 0; i < ndof; ++i)
              {
                dshape (i * 4 + d, imip) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (i * 4 + d, imip)
                      += (localmat)[2][j] * pol[(localmat)[1][j]] * scale[d];
              }
          }
      }
  }

  template <>
  void
  BlockMappedElement<1>::CalcMappedDDShape (const BaseMappedIntegrationPoint &,
                                            BareSliceMatrix<>) const
  {
    ;
  }

  template <>
  void BlockMappedElement<2>::CalcMappedDDShape (
      const BaseMappedIntegrationPoint &bmip, BareSliceMatrix<> hddshape) const
  {
    auto ddshape = hddshape.AddSize (this->ndof, 2 * 2);

    // auto & mip = static_cast<const MappedIntegrationPoint<2,2> &> (bmip);
    Vec<2> cpoint = bmip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 2 * (order + 1) + 1);
    mem2[0] = 0;
    double *polxt[2];
    for (int d = 0; d < 2; d++)
      {
        polxt[d] = &mem2[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d1 = 0; d1 < 2; d1++)
      {
        CSR localmat = localmats[d1];
        for (int d2 = 0; d2 < 2; d2++)
          {
            Vector<double> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                pol[ii++] = (d2 == 0 ? i : j) * polxt[0][i - (d2 == 0)]
                            * polxt[1][j - (d2 == 1)];

            for (int i = 0; i < this->ndof; ++i)
              {
                ddshape (i, d2 * 2 + d1) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  ddshape (i, d2 * 2 + d1) += (localmat)[2][j]
                                              * pol[(localmat)[1][j]]
                                              * scale[d1] * scale[d2];
              }
          }
      }
  }

  template <>
  void BlockMappedElement<3>::CalcMappedDDShape (
      const BaseMappedIntegrationPoint &bmip, BareSliceMatrix<> hddshape) const
  {
    auto ddshape = hddshape.AddSize (this->ndof, 3 * 3);
    // auto & mip = static_cast<const MappedIntegrationPoint<2,2> &> (bmip);
    Vec<3> cpoint = bmip.GetPoint ();
    cpoint -= shift;
    vtimes (cpoint, scale);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 3 * (order + 1) + 1);
    mem2[0] = 0;
    double *polxt[3];
    for (int d = 0; d < 3; d++)
      {
        polxt[d] = &mem2[d * (order + 1) + 1];
        Monomial (order, cpoint[d], polxt[d]);
      }

    for (int d1 = 0; d1 < 3; d1++)
      {
        CSR localmat = localmats[d1];
        for (int d2 = 0; d2 < 3; d2++)
          {
            Vector<double> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                for (int k = 0; k <= order - i - j; k++)
                  pol[ii++] = (d2 == 0 ? i : (d2 == 1 ? j : k))
                              * polxt[0][i - (d2 == 0)]
                              * polxt[1][j - (d2 == 1)]
                              * polxt[2][k - (d2 == 2)];

            for (int i = 0; i < this->ndof; ++i)
              {
                ddshape (i, d2 * 3 + d1) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  ddshape (i, d2 * 3 + d1) += (localmat)[2][j]
                                              * pol[(localmat)[1][j]]
                                              * scale[d1] * scale[d2];
              }
          }
      }
  }

  template <>
  void
  BlockMappedElement<4>::CalcMappedDDShape (const BaseMappedIntegrationPoint &,
                                            BareSliceMatrix<>) const
  {
    ;
  }

  template class BlockMappedElement<1>;
  template class BlockMappedElement<2>;
  template class BlockMappedElement<3>;
  template class BlockMappedElement<4>;

}
