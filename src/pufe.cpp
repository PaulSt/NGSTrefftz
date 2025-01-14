#include "pufe.hpp"

namespace ngfem
{

  template <int D> string PUFElement<D>::ClassName () const
  {
    return "PUFElement";
  }

  template <int D>
  void PUFElement<D>::CalcDShape (const BaseMappedIntegrationRule &mir,
                                  BareSliceMatrix<> dshapes) const
  {
    for (size_t i = 0; i < mir.Size (); i++)
      CalcDShape (mir[i], dshapes.Cols (i * D, (i + 1) * D));
  }

  template <int D>
  void PUFElement<D>::CalcMappedDShape (const BaseMappedIntegrationRule &mir,
                                        BareSliceMatrix<> dshapes) const
  {
    for (size_t i = 0; i < mir.Size (); i++)
      CalcMappedDShape (mir[i], dshapes.Cols (i * D, (i + 1) * D));
  }

  template <int D>
  Vec<D> PUFElement<D>::EvaluateGrad (const BaseMappedIntegrationPoint &ip,
                                      BareSliceVector<double> x) const
  {
    MatrixFixWidth<D> dshape (ndof);
    CalcDShape (ip, dshape);
    Vec<D> grad = Trans (dshape) * x;
    return grad;
  }

  template <int D>
  void PUFElement<D>::EvaluateGrad (const BaseMappedIntegrationRule &ir,
                                    BareSliceVector<double> coefs,
                                    FlatMatrixFixWidth<D, double> vals) const
  {
    for (size_t i = 0; i < ir.Size (); i++)
      vals.Row (i) = EvaluateGrad (ir[i], coefs);
  }

  template <int D>
  void PUFElement<D>::EvaluateGradTrans (const BaseMappedIntegrationRule &ir,
                                         FlatMatrixFixWidth<D, double> vals,
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
  void PUFElement<D>::EvaluateGradTrans (const BaseMappedIntegrationRule &,
                                         SliceMatrix<>, SliceMatrix<>) const
  {
    cout << "EvalGradTrans not overloaded" << endl;
  }

  ////////////////////////////////////////////////////////
  /////////////// CalcShape implementation ///////////////
  ////////////////////////////////////////////////////////

  template <>
  void PUFElement<1>::CalcShape (const SIMD_BaseMappedIntegrationRule &,
                                 BareSliceMatrix<SIMD<double>>) const
  {
    throw ExceptionNOSIMD ("nosimd");
  }

  template <>
  void PUFElement<2>::CalcShape (const SIMD_BaseMappedIntegrationRule &,
                                 BareSliceMatrix<SIMD<double>>) const
  {
    throw ExceptionNOSIMD ("nosimd");
  }

  template <>
  void PUFElement<3>::CalcShape (const SIMD_BaseMappedIntegrationRule &,
                                 BareSliceMatrix<SIMD<double>>) const
  {
    throw ExceptionNOSIMD ("nosimd");
  }

  template <>
  void PUFElement<1>::CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                                  BareSliceMatrix<SIMD<double>>) const
  {
    throw ExceptionNOSIMD ("nosimd");
  }

  template <>
  void PUFElement<2>::CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                                  BareSliceMatrix<SIMD<double>>) const
  {
    throw ExceptionNOSIMD ("nosimd");
  }

  template <>
  void PUFElement<3>::CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                                  BareSliceMatrix<SIMD<double>>) const
  {
    throw ExceptionNOSIMD ("nosimd");
  }

  /////////////// non-simd

  template <>
  void PUFElement<1>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                 BareSliceVector<> shape) const
  {
    auto ip = mip.IP ();
    double lam[2] = { ip (0), 1 - ip (0) };
    for (int v = 0; v < 2; v++)
      {
        Vec<1> cpoint = mip.GetPoint ();
        cpoint -= elvertices[v];
        cpoint *= (1.0 / elsizes[v]);
        //  calc 1 dimensional monomial basis
        STACK_ARRAY (double, mem, (order + 1));
        // double *polxt[1];
        // polxt[0] = &mem[0];
        FlatVector<double> polxt (order + 1, &mem[0]);
        Monomial (order, cpoint[0], &mem[0]);
        // TB*monomials for trefftz shape fcts
        int nlbasis = localmat[0].Size () - 1;

        for (int i = 0; i < nlbasis; ++i)
          {
            shape (nlbasis * v + i) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (nlbasis * v + i)
                  += (localmat)[2][j] * polxt[(localmat)[1][j]];
            shape (nlbasis * v + i) *= lam[v];
          }
      }
  }

  template <>
  void PUFElement<2>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                 BareSliceVector<> shape) const
  {
    auto ip = mip.IP ();
    double lam[3] = { ip (0), ip (1), 1 - ip (0) - ip (1) };
    for (int v = 0; v < 3; v++)
      {
        // auto & smir = static_cast<const SIMD_BaseMappedIntegrationRuleD>&>
        // (mir);
        Vec<2> cpoint = mip.GetPoint ();
        cpoint -= elvertices[v];
        cpoint *= (1.0 / elsizes[v]);
        //  calc 1 dimensional monomial basis
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
        int nlbasis = localmat[0].Size () - 1;

        // cout << "npoly: " << npoly << endl;
        // cout << "nlbasis: " << nlbasis << endl;
        //  cout << "npoly: " << npoly << endl;
        //  cout << "order: " << order << endl;
        //  cout << "ndof: " << this->ndof << endl;
        //  cout << nlbasis * 3 << endl;
        //  cout << "npoly: " << npoly << endl;
        //  cout << "localmat[0].Size(): " << localmat[0].Size() << endl;
        //  cout << "elvertex: " << v << " : " << elvertices[v] << endl;
        //
        for (int i = 0; i < nlbasis; ++i)
          {
            shape (nlbasis * v + i) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (nlbasis * v + i)
                  += (localmat)[2][j] * pol[(localmat)[1][j]];
            shape (nlbasis * v + i) *= lam[v];
          }
      }
  }

  template <>
  void PUFElement<3>::CalcShape (const BaseMappedIntegrationPoint &,
                                 BareSliceVector<>) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void PUFElement<1>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                  BareSliceMatrix<> dshape) const
  {
    auto ip = mip.IP ();
    double lam[2] = { ip (0), 1 - ip (0) };
    double dlam[4] = { 1, -1 };
    FlatMatrix<> dlam_mat (2, 1, dlam);
    for (int v = 0; v < 2; v++)
      {
        Vec<1> cpoint = mip.GetPoint ();
        cpoint -= elvertices[v];
        cpoint *= (1.0 / elsizes[v]);
        // cpoint[1] *= c;
        //  +1 size to avoid undefined behavior taking deriv, getting [-1]
        //  entry

        auto JI = static_cast<const MappedIntegrationPoint<1, 1> &> (mip)
                      .GetJacobianInverse ();
        auto dlam_mapped = dlam_mat * JI;

        STACK_ARRAY (double, mem, (order + 1) + 1);
        mem[0] = 0;
        FlatVector<double> polxt (order + 2, &mem[0]);
        Monomial (order, cpoint[0], &mem[1]);

        int nlbasis = localmat[0].Size () - 1;

        Vector<double> shape (2 * nlbasis);
        // TB*monomials for trefftz shape fcts
        for (int i = 0; i < nlbasis; ++i)
          {
            shape (nlbasis * v + i) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (nlbasis * v + i)
                  += (localmat)[2][j] * polxt[(localmat)[1][j] + 1];
            // shape (nlbasis * v + i) *= lam[v];
          }

        for (int i = 0; i < nlbasis; ++i)
          {
            dshape (nlbasis * v + i, 0) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (nlbasis * v + i, 0)
                  += (localmat)[2][j] * polxt[(localmat)[1][j]]
                     * (localmat)[1][j] * (1.0 / elsizes[v]);
            dshape (nlbasis * v + i, 0) *= lam[v];

            dshape (nlbasis * v + i, 0)
                += dlam_mapped (v, 0) * shape (nlbasis * v + i);
            // static_pointer_cast<const MappedIntegrationPoint<2,2>>(mip);
            // //.GetJacobianInverse(); mip.GetJacobianInverse()
          }
      }
  }

  template <>
  void PUFElement<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                  BareSliceMatrix<> dshape) const
  {
    auto ip = mip.IP ();
    double lam[3] = { ip (0), ip (1), 1 - ip (0) - ip (1) };
    // double dlam[3][2] = { { 1, 0 }, { 0, 1 }, { -1, -1 } };
    double dlam[6] = { 1, 0, 0, 1, -1, -1 };
    FlatMatrix<> dlam_mat (3, 2, dlam);
    for (int v = 0; v < 3; v++)
      {
        Vec<2> cpoint = mip.GetPoint ();
        cpoint -= elvertices[v];
        cpoint *= (1.0 / elsizes[v]);
        // cpoint[1] *= c;
        //  +1 size to avoid undefined behavior taking deriv, getting [-1]
        //  entry

        auto JI = static_cast<const MappedIntegrationPoint<2, 2> &> (mip)
                      .GetJacobianInverse ();
        auto dlam_mapped = dlam_mat * JI;

        STACK_ARRAY (double, mem, 2 * (order + 1) + 1);
        mem[0] = 0;
        double *polxt[2];
        for (int d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (order + 1) + 1];
            Monomial (order, cpoint[d], polxt[d]);
          }

        int nlbasis = localmat[0].Size () - 1;

        Vector<double> pol1 (npoly);
        Vector<double> shape (3 * nlbasis);
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            pol1[ii++] = polxt[0][i] * polxt[1][j];
        // TB*monomials for trefftz shape fcts
        for (int i = 0; i < nlbasis; ++i)
          {
            shape (nlbasis * v + i) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              shape (nlbasis * v + i)
                  += (localmat)[2][j] * pol1[(localmat)[1][j]];
            // shape (nlbasis * v + i) *= lam[v];
          }

        for (int d = 0; d < 2; d++)
          {
            Vector<double> pol (npoly);
            for (int i = 0, ii = 0; i <= order; i++)
              for (int j = 0; j <= order - i; j++)
                pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                            * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

            for (int i = 0; i < nlbasis; ++i)
              {
                dshape (nlbasis * v + i, d) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  dshape (nlbasis * v + i, d) += (localmat)[2][j]
                                                 * pol[(localmat)[1][j]]
                                                 * (1.0 / elsizes[v]);
                dshape (nlbasis * v + i, d) *= lam[v];

                dshape (nlbasis * v + i, d)
                    += dlam_mapped (v, d) * shape (nlbasis * v + i);
                // static_pointer_cast<const MappedIntegrationPoint<2,2>>(mip);
                // //.GetJacobianInverse(); mip.GetJacobianInverse()
              }
          }
      }
  }

  template <>
  void PUFElement<3>::CalcDShape (const BaseMappedIntegrationPoint &,
                                  BareSliceMatrix<>) const
  {
    cout << "dim not implemented" << endl;
  }

  template class PUFElement<1>;
  template class PUFElement<2>;
  template class PUFElement<3>;

}
