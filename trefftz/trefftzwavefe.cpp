#include "trefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
  template <int D>
  TrefftzWaveFE<D>::TrefftzWaveFE (int aord, float ac, Vec<D> aelcenter,
                                   double aelsize, ELEMENT_TYPE aeltype)
      : ScalarMappedElement<D> (BinCoeff (D - 1 + aord, aord)
                                    + BinCoeff (D - 1 + aord - 1, aord - 1),
                                aord),
        ord (aord), c (ac), nbasis (BinCoeff (D - 1 + ord, ord)
                                    + BinCoeff (D - 1 + ord - 1, ord - 1)),
        npoly (BinCoeff (D + ord, ord)), elcenter (aelcenter),
        elsize (aelsize), eltype (aeltype)
  {
    ;
  }

  template <>
  void TrefftzWaveFE<1>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzWaveFE<2>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[1] *= c;
        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 2 * (ord + 1));
        Vec<2, SIMD<double> *> polxt;
        for (size_t d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (ord + 1)];
            Monomial (ord, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            pol[ii++] = polxt[0][i] * polxt[1][j];
        // TB*monomials for trefftz shape fcts
        const CSR *localmat = TrefftzWaveBasis<2>::getInstance ().TB (ord);
        for (int i = 0; i < nbasis; ++i)
          {
            shape (i, imip) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              shape (i, imip) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
  }

  template <>
  void TrefftzWaveFE<3>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[2] *= c;
        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 3 * (ord + 1));
        Vec<3, SIMD<double> *> polxt;
        for (size_t d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (ord + 1)];
            Monomial (ord, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            for (size_t k = 0; k <= ord - i - j; k++)
              pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];
        // TB*monomials for trefftz shape fcts
        const CSR *localmat = TrefftzWaveBasis<3>::getInstance ().TB (ord);
        for (int i = 0; i < nbasis; ++i)
          {
            shape (i, imip) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              shape (i, imip) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
  }

  template <>
  void TrefftzWaveFE<4>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    // auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&>
    // (mir);
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<4, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[3] *= c;
        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 4 * (ord + 1));
        Vec<4, SIMD<double> *> polxt;
        for (size_t d = 0; d < 4; d++)
          {
            polxt[d] = &mem[d * (ord + 1)];
            Monomial (ord, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            for (size_t k = 0; k <= ord - i - j; k++)
              for (size_t l = 0; l <= ord - i - j - k; l++)
                pol[ii++]
                    = polxt[0][i] * polxt[1][j] * polxt[2][k] * polxt[3][l];
        // TB*monomials for trefftz shape fcts
        const CSR *localmat = TrefftzWaveBasis<4>::getInstance ().TB (ord);
        for (int i = 0; i < nbasis; ++i)
          {
            shape (i, imip) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              shape (i, imip) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
  }

  template <>
  void
  TrefftzWaveFE<1>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
  }

  template <>
  void
  TrefftzWaveFE<2>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[1] *= c;

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 2 * (ord + 1) + 1);
        mem[0] = 0;
        Vec<2, SIMD<double> *> polxt;
        for (size_t d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (ord + 1) + 1];
            Monomial (ord, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 2; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
              for (size_t j = 0; j <= ord - i; j++)
                pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                            * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

            const CSR *localmat = TrefftzWaveBasis<2>::getInstance ().TB (ord);
            for (int i = 0; i < nbasis; ++i)
              {
                dshape (i * 2 + d, imip) = 0.0;
                for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
                  dshape (i * 2 + d, imip)
                      += (*localmat)[2][j] * pol[(*localmat)[1][j]]
                         * (d == 1 ? c : 1) * (2.0 / elsize);
              }
          }
      }
    dshape *= (2.0 / elsize); // inner derivative
  }

  template <>
  void
  TrefftzWaveFE<3>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[2] *= c;

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 3 * (ord + 1) + 1);
        mem[0] = 0;
        Vec<3, SIMD<double> *> polxt;
        for (size_t d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (ord + 1) + 1];
            Monomial (ord, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 3; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
              for (size_t j = 0; j <= ord - i; j++)
                for (size_t k = 0; k <= ord - i - j; k++)
                  pol[ii++] = (d == 0 ? i : (d == 1 ? j : (d == 2 ? k : 0)))
                              * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                              * polxt[2][k - (d == 2)];

            const CSR *localmat = TrefftzWaveBasis<3>::getInstance ().TB (ord);
            for (int i = 0; i < nbasis; ++i)
              {
                dshape (i * 3 + d, imip) = 0.0;
                for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
                  dshape (i * 3 + d, imip)
                      += (*localmat)[2][j] * pol[(*localmat)[1][j]]
                         * (d == 2 ? c : 1) * (2.0 / elsize);
              }
          }
      }
    dshape *= (2.0 / elsize); // inner derivative
  }

  template <>
  void
  TrefftzWaveFE<4>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<4, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[3] *= c;

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 4 * (ord + 1) + 1);
        mem[0] = 0;
        Vec<4, SIMD<double> *> polxt;
        for (size_t d = 0; d < 4; d++)
          {
            polxt[d] = &mem[d * (ord + 1) + 1];
            Monomial (ord, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 4; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
              for (size_t j = 0; j <= ord - i; j++)
                for (size_t k = 0; k <= ord - i - j; k++)
                  for (size_t l = 0; l <= ord - i - j - k; l++)
                    pol[ii++]
                        = (d == 0 ? i
                                  : (d == 1 ? j
                                            : (d == 2 ? k : (d == 3 ? l : 0))))
                          * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                          * polxt[2][k - (d == 2)] * polxt[3][l - (d == 3)];

            const CSR *localmat = TrefftzWaveBasis<4>::getInstance ().TB (ord);
            for (int i = 0; i < nbasis; ++i)
              {
                dshape (i * 4 + d, imip) = 0.0;
                for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
                  dshape (i * 4 + d, imip)
                      += (*localmat)[2][j] * pol[(*localmat)[1][j]]
                         * (d == 3 ? c : 1) * (2.0 / elsize);
              }
          }
      }
    // dshape.AddSize(nbasis*4,smir.Size()) *= (2.0/elsize); //inner derivative
  }

  /////////////// non-simd

  template <>
  void TrefftzWaveFE<1>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzWaveFE<2>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    // auto & smir = static_cast<const SIMD_BaseMappedIntegrationRuleD>&>
    // (mir);
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[1] *= c;
    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 2 * (ord + 1));
    double *polxt[2];
    for (size_t d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (ord + 1)];
        Monomial (ord, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (size_t i = 0, ii = 0; i <= ord; i++)
      for (size_t j = 0; j <= ord - i; j++)
        pol[ii++] = polxt[0][i] * polxt[1][j];
    // TB*monomials for trefftz shape fcts
    const CSR *localmat = TrefftzWaveBasis<2>::getInstance ().TB (ord);
    for (int i = 0; i < nbasis; ++i)
      {
        shape (i) = 0.0;
        for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
          shape (i) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
      }
  }

  template <>
  void TrefftzWaveFE<3>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    Vec<3> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[2] *= c;
    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 3 * (ord + 1));
    // double* polxt[3];
    Vec<3, double *> polxt;
    for (size_t d = 0; d < 3; d++)
      {
        polxt[d] = &mem[d * (ord + 1)];
        Monomial (ord, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (size_t i = 0, ii = 0; i <= ord; i++)
      for (size_t j = 0; j <= ord - i; j++)
        for (size_t k = 0; k <= ord - i - j; k++)
          pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];
    // TB*monomials for trefftz shape fcts
    const CSR *localmat = TrefftzWaveBasis<3>::getInstance ().TB (ord);
    for (int i = 0; i < nbasis; ++i)
      {
        shape (i) = 0.0;
        for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
          shape (i) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
      }
  }

  template <>
  void TrefftzWaveFE<4>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    Vec<4> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[3] *= c;
    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 4 * (ord + 1));
    double *polxt[4];
    for (size_t d = 0; d < 4; d++)
      {
        polxt[d] = &mem[d * (ord + 1)];
        Monomial (ord, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (size_t i = 0, ii = 0; i <= ord; i++)
      for (size_t j = 0; j <= ord - i; j++)
        for (size_t k = 0; k <= ord - i - j; k++)
          for (size_t l = 0; l <= ord - i - j - k; l++)
            pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k] * polxt[3][l];
    // TB*monomials for trefftz shape fcts
    const CSR *localmat = TrefftzWaveBasis<4>::getInstance ().TB (ord);
    for (int i = 0; i < nbasis; ++i)
      {
        shape (i) = 0.0;
        for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
          shape (i) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
      }
  }

  template <>
  void TrefftzWaveFE<1>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzWaveFE<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[1] *= c;
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 2 * (ord + 1) + 1);
    mem[0] = 0;
    double *polxt[2];
    for (size_t d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (ord + 1) + 1];
        Monomial (ord, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 2; d++)
      {
        Vector<double> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                        * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

        const CSR *localmat = TrefftzWaveBasis<2>::getInstance ().TB (ord);
        for (int i = 0; i < nbasis; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              dshape (i, d) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
    dshape *= (2.0 / elsize); // inner derivative
    dshape.Col (1) *= c;      // inner derivative
  }

  template <>
  void TrefftzWaveFE<3>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    Vec<3> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[2] *= c;
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 3 * (ord + 1) + 1);
    mem[0] = 0;
    // double* polxt[4];
    Vec<3, double *> polxt;
    for (size_t d = 0; d < 3; d++)
      {
        polxt[d] = &mem[d * (ord + 1) + 1];
        Monomial (ord, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 3; d++)
      {
        Vector<double> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            for (size_t k = 0; k <= ord - i - j; k++)
              pol[ii++] = (d == 0 ? i : (d == 1 ? j : (d == 2 ? k : 0)))
                          * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                          * polxt[2][k - (d == 2)];

        const CSR *localmat = TrefftzWaveBasis<3>::getInstance ().TB (ord);
        for (int i = 0; i < nbasis; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              dshape (i, d) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
    dshape *= (2.0 / elsize); // inner derivative
    dshape.Col (2) *= c;      // inner derivative
  }

  template <>
  void TrefftzWaveFE<4>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    Vec<4> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    cpoint[3] *= c;
    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 4 * (ord + 1) + 1);
    mem[0] = 0;
    double *polxt[4];
    for (size_t d = 0; d < 4; d++)
      {
        polxt[d] = &mem[d * (ord + 1) + 1];
        Monomial (ord, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 4; d++)
      {
        Vector<double> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            for (size_t k = 0; k <= ord - i - j; k++)
              for (size_t l = 0; l <= ord - i - j - k; l++)
                pol[ii++]
                    = (d == 0 ? i
                              : (d == 1 ? j : (d == 2 ? k : (d == 3 ? l : 0))))
                      * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)]
                      * polxt[2][k - (d == 2)] * polxt[3][l - (d == 3)];

        const CSR *localmat = TrefftzWaveBasis<4>::getInstance ().TB (ord);
        for (int i = 0; i < nbasis; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              dshape (i, d) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
    dshape *= (2.0 / elsize); // inner derivative
    dshape.Col (3) *= c;      // inner derivative
  }

  template class TrefftzWaveFE<1>;
  template class TrefftzWaveFE<2>;
  template class TrefftzWaveFE<3>;
  template class TrefftzWaveFE<4>;

  template <int D> const CSR *TrefftzWaveBasis<D>::TB (int ord)
  {
    const CSR *tb = &tbstore[ord];
    return tb;
  }

  template <int D> void TrefftzWaveBasis<D>::CreateTB (int ord)
  {
    cout << "creating tp store for order " << ord << endl;
    // if (tbstore.Size() < ord)
    //{
    int oldsize = tbstore.Size ();
    tbstore.SetSize (ord + 1);
    for (int i = oldsize; i <= ord; i++)
      tbstore[i] = CSR ();

    // if ( tbstore[ord][0].Size() == 0)
    //{
    const int nbasis
        = (BinCoeff (D - 1 + ord, ord) + BinCoeff (D - 1 + ord - 1, ord - 1));
    const int npoly = (BinCoeff (D + ord, ord));
    Matrix<> trefftzbasis (nbasis, npoly);
    trefftzbasis = 0;
    Vec<D, int> coeff = 0;
    int count = 0;
    for (int b = 0; b < nbasis; b++)
      {
        int tracker = 0;
        TB_inner (ord, trefftzbasis, coeff, b, D, tracker);
      }
    MatToCSR (trefftzbasis, tbstore[ord]);
    //}
    //}

    // if ( tbstore[ord][0].Size() == 0)
    if (tbstore.Size () < ord)
      {
        stringstream str;
        str << "failed to generate trefftz basis of order " << ord << endl;
        throw Exception (str.str ());
      }
  }

  template <int D>
  void TrefftzWaveBasis<D>::TB_inner (int ord, Matrix<> &trefftzbasis,
                                      Vec<D, int> coeffnum, int basis, int dim,
                                      int &tracker)
  {
    if (dim > 0)
      {
        while (coeffnum (dim - 1) <= ord)
          {
            TB_inner (ord, trefftzbasis, coeffnum, basis, dim - 1, tracker);
            coeffnum (dim - 1)++;
          }
      }
    else
      {
        int sum = 0;
        for (int i = 0; i < D; i++)
          sum += coeffnum (i);
        if (sum <= ord)
          {
            if (tracker >= 0)
              tracker++;
            int indexmap = IndexMap2 (coeffnum, ord);
            if ((coeffnum (D - 1) == 0 || coeffnum (D - 1) == 1)
                && tracker > basis)
              {
                trefftzbasis (basis, indexmap) = 1;
                tracker = -1;
              }
            else if (coeffnum (D - 1) > 1)
              {
                int k = coeffnum (D - 1);
                for (int m = 0; m < D - 1; m++) // rekursive sum
                  {
                    Vec<D, int> get_coeff = coeffnum;
                    get_coeff[D - 1] = get_coeff[D - 1] - 2;
                    get_coeff[m] = get_coeff[m] + 2;
                    trefftzbasis (basis, indexmap)
                        += (coeffnum (m) + 1) * (coeffnum (m) + 2)
                           * trefftzbasis (basis, IndexMap2 (get_coeff, ord));
                  }
                trefftzbasis (basis, indexmap) *= 1.0 / (k * (k - 1));
              }
          }
      }
  }

  template <int D>
  int TrefftzWaveBasis<D>::IndexMap2 (Vec<D, int> index, int ord)
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

  template class TrefftzWaveBasis<1>;
  template class TrefftzWaveBasis<2>;
  template class TrefftzWaveBasis<3>;
  template class TrefftzWaveBasis<4>;

  void MatToCSR (Matrix<> mat, CSR &sparsemat)
  {
    int spsize = 0;
    for (int i = 0; i < mat.Height (); i++)
      {
        // guarantee sparsemat[0].Size() to be nbasis,
        // could be shorter but easier for mat*vec
        sparsemat[0].Append (spsize);
        for (int j = 0; j < mat.Width (); j++)
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NGS_PYTHON
void ExportTrefftzElement (py::module m)
{
  // py::class_<TrefftzWaveFE<3>, shared_ptr<TrefftzWaveFE<3>>, FiniteElement>
  // 	(m, "TrefftzWaveFE3", "Trefftz space for wave eq")
  // 	.def(py::init<>())
  // 	;
  // py::class_<TrefftzWaveFE<2>, shared_ptr<TrefftzWaveFE<2>>, FiniteElement>
  // 	(m, "TrefftzWaveFE2", "Trefftz space for wave eq")
  // 	.def(py::init<>())
  // 	;
  // py::class_<TrefftzWaveFE<1>, shared_ptr<TrefftzWaveFE<1>>, FiniteElement>
  // 	(m, "TrefftzWaveFE1", "Trefftz space for wave eq")
  // 	.def(py::init<>())
  // 	;
}
#endif // NGS_PYTHON
