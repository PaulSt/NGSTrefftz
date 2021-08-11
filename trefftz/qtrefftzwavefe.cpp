#include "qtrefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"
#include "trefftzwavefe.hpp"

#include <ctime>

namespace ngfem
{

  template <>
  void QTrefftzWaveFE<1>::CalcDDSpecialShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape,
      BareSliceMatrix<SIMD<double>> wavespeed,
      BareSliceMatrix<SIMD<double>> mu) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (1.0 / elsize);

        STACK_ARRAY (SIMD<double>, mem, 2 * (ord + 1) + 2);
        mem[0] = 0;
        mem[1] = 0;
        Vec<2, SIMD<double> *> polxt;
        for (size_t d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (ord + 1) + 2];
            Monomial (ord, cpoint[d], polxt[d]);
          }

        Vector<SIMD<double>> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            pol[ii++] = (i * (i - 1) * polxt[0][i - 2] * polxt[1][j]
                         - j * (j - 1) * polxt[0][i] * polxt[1][j - 2]
                               * wavespeed (0, imip))
                        * mu (0, imip);

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i * 2, imip) = 0.0;
            dshape (i * 2 + 1, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i * 2 + 1, imip) += (localmat)[2][j]
                                          * pol[(localmat)[1][j]]
                                          * pow (1.0 / elsize, 2);
          }
      }
  }

  template <>
  void QTrefftzWaveFE<2>::CalcDDSpecialShape (
      const SIMD_BaseMappedIntegrationRule &smir,
      BareSliceMatrix<SIMD<double>> dshape,
      BareSliceMatrix<SIMD<double>> wavespeed,
      BareSliceMatrix<SIMD<double>> mu) const
  {
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<3, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (1.0 / elsize);

        STACK_ARRAY (SIMD<double>, mem, 3 * (ord + 1) + 2);
        mem[0] = 0;
        mem[1] = 0;
        Vec<3, SIMD<double> *> polxt;
        for (size_t d = 0; d < 3; d++)
          {
            polxt[d] = &mem[d * (ord + 1) + 2];
            Monomial (ord, cpoint[d], polxt[d]);
          }

        Vector<SIMD<double>> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            for (size_t k = 0; k <= ord - i - j; k++)
              pol[ii++]
                  = (i * (i - 1) * polxt[0][i - 2] * polxt[1][j] * polxt[2][k]
                     + j * (j - 1) * polxt[0][i] * polxt[1][j - 2]
                           * polxt[2][k]
                     - k * (k - 1) * polxt[0][i] * polxt[1][j]
                           * polxt[2][k - 2] * wavespeed (0, imip))
                    * mu (0, imip);

        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i * 3, imip) = 0.0;
            dshape (i * 3 + 1, imip) = 0.0;
            dshape (i * 3 + 2, imip) = 0.0;
            for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
              dshape (i * 3 + 2, imip) += (localmat)[2][j]
                                          * pol[(localmat)[1][j]]
                                          * pow (1.0 / elsize, 2);
          }
      }
  }

  template <>
  void
  QTrefftzWaveFE<1>::CalcMappedDDShape (const BaseMappedIntegrationPoint &bmip,
                                        BareSliceMatrix<> hddshape) const
  {
    auto ddshape = hddshape.AddSize (this->ndof, 2 * 2);

    // auto & mip = static_cast<const MappedIntegrationPoint<2,2> &> (bmip);
    Vec<2> cpoint = bmip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (1.0 / elsize);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 2 * (ord + 1) + 1);
    mem2[0] = 0;
    double *polxt[2];
    for (size_t d = 0; d < 2; d++)
      {
        polxt[d] = &mem2[d * (ord + 1) + 1];
        Monomial (ord, cpoint[d], polxt[d]);
      }

    for (int d1 = 0; d1 < 2; d1++)
      {
        for (int d2 = 0; d2 < 2; d2++)
          {
            Vector<double> pol (npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
              for (size_t j = 0; j <= ord - i; j++)
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
                                              * pow (1.0 / elsize, 2);
              }
          }
      }
  }

  template <>
  void
  QTrefftzWaveFE<2>::CalcMappedDDShape (const BaseMappedIntegrationPoint &bmip,
                                        BareSliceMatrix<> hddshape) const
  {
    auto ddshape = hddshape.AddSize (this->ndof, 3 * 3);
    // auto & mip = static_cast<const MappedIntegrationPoint<2,2> &> (bmip);
    Vec<3> cpoint = bmip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (1.0 / elsize);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 3 * (ord + 1) + 1);
    mem2[0] = 0;
    double *polxt[3];
    for (size_t d = 0; d < 3; d++)
      {
        polxt[d] = &mem2[d * (ord + 1) + 1];
        Monomial (ord, cpoint[d], polxt[d]);
      }

    for (int d1 = 0; d1 < 3; d1++)
      {
        for (int d2 = 0; d2 < 3; d2++)
          {
            Vector<double> pol (npoly);
            for (int i = 0, ii = 0; i <= ord; i++)
              for (int j = 0; j <= ord - i; j++)
                for (int k = 0; k <= ord - i - j; k++)
                  pol[ii++] = (d1 == d2 ? (d1 == 0 ? i * (i - 1)
                                                   : (d1 == 1 ? j * (j - 1)
                                                              : k * (k - 1)))
                                        : (int[3]){ i, j, k }[d1]
                                              * (int[3]){ i, j, k }[d2])
                              * polxt[0][i - (d1 == 0) - (d2 == 0)]
                              * polxt[1][j - (d1 == 1) - (d2 == 1)]
                              * polxt[2][k - (d1 == 2) - (d2 == 2)];

            for (int i = 0; i < this->ndof; ++i)
              {
                ddshape (i, d2 * 3 + d1) = 0.0;
                for (int j = (localmat)[0][i]; j < (localmat)[0][i + 1]; ++j)
                  ddshape (i, d2 * 3 + d1) += (localmat)[2][j]
                                              * pol[(localmat)[1][j]]
                                              * pow (1.0 / elsize, 2);
              }
          }
      }
  }

  template class QTrefftzWaveFE<1>;
  template class QTrefftzWaveFE<2>;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D>
  CSR QTrefftzWaveBasis<D>::TB (int ord, Vec<D + 1> ElCenter,
                                Matrix<shared_ptr<CoefficientFunction>> GGder,
                                Matrix<shared_ptr<CoefficientFunction>> BBder,
                                double elsize, int basistype)
  {
    lock_guard<mutex> lock (gentrefftzbasis);
    string encode = to_string (ord);
    for (int i = 0; i < D; i++)
      encode += to_string (ElCenter[i]);

    if (gtbstore[encode][0].Size () == 0)
      {
        IntegrationRule ir (D == 3 ? ET_TET : D == 2 ? ET_TRIG : ET_SEGM, 0);
        Mat<D, D> dummy;
        FE_ElementTransformation<D, D> et (D == 3   ? ET_TET
                                           : D == 2 ? ET_TRIG
                                                    : ET_SEGM,
                                           dummy);
        MappedIntegrationPoint<D, D> mip (ir[0], et, 0);
        for (int i = 0; i < D; i++)
          mip.Point ()[i] = ElCenter[i];

        Matrix<> BB (ord, (ord - 1) * (D == 2) + 1);
        Matrix<> GG (ord - 1, (ord - 2) * (D == 2) + 1);
        for (int ny = 0; ny <= (ord - 1) * (D == 2); ny++)
          {
            for (int nx = 0; nx <= ord - 1; nx++)
              {
                double fac = (factorial (nx) * factorial (ny));
                BB (nx, ny) = BBder (nx, ny)->Evaluate (mip) / fac
                              * pow (elsize, nx + ny);
                if (nx < ord - 1 && ny < ord - 1)
                  GG (nx, ny) = GGder (nx, ny)->Evaluate (mip) / fac
                                * pow (elsize, nx + ny);
              }
          }

        const int nbasis
            = (BinCoeff (D + ord, ord) + BinCoeff (D + ord - 1, ord - 1));
        const int npoly = BinCoeff (D + 1 + ord, ord);
        Matrix<> qbasis (nbasis, npoly);
        qbasis = 0;

        for (int t = 0, basisn = 0; t < 2; t++)
          for (int x = 0; x <= ord - t; x++)
            for (int y = 0; y <= (ord - x - t) * (D == 2); y++)
              {
                Vec<D + 1, int> index;
                index[D] = t;
                index[0] = x;
                if (D == 2)
                  index[1] = y;
                qbasis (basisn++, TrefftzWaveBasis<D>::IndexMap2 (index, ord))
                    = 1.0;
              }

        for (int basisn = 0; basisn < nbasis; basisn++)
          {
            for (int ell = 0; ell < ord - 1; ell++)
              {
                for (int t = 0; t <= ell; t++)
                  {
                    for (int x = (D == 1 ? ell - t : 0); x <= ell - t; x++)
                      {
                        int y = ell - t - x;
                        Vec<D + 1, int> index;
                        index[1] = y;
                        index[0] = x;
                        index[D] = t + 2;
                        double *newcoeff = &qbasis (
                            basisn,
                            TrefftzWaveBasis<D>::IndexMap2 (index, ord));
                        *newcoeff = 0;

                        for (int betax = 0; betax <= x; betax++)
                          for (int betay = (D == 2) ? 0 : y; betay <= y;
                               betay++)
                            {
                              index[1] = betay;
                              index[0] = betax + 1;
                              index[D] = t;
                              int getcoeffx = TrefftzWaveBasis<D>::IndexMap2 (
                                  index, ord);
                              index[1] = betay + 1;
                              index[0] = betax;
                              index[D] = t;
                              int getcoeffy = TrefftzWaveBasis<D>::IndexMap2 (
                                  index, ord);
                              index[1] = betay;
                              index[0] = betax + 2;
                              index[D] = t;
                              int getcoeffxx = TrefftzWaveBasis<D>::IndexMap2 (
                                  index, ord);
                              index[1] = betay + 2;
                              index[0] = betax;
                              index[D] = t;
                              int getcoeffyy = TrefftzWaveBasis<D>::IndexMap2 (
                                  index, ord);

                              *newcoeff
                                  += (betax + 2) * (betax + 1)
                                         / ((t + 2) * (t + 1) * GG (0))
                                         * BB (x - betax, y - betay)
                                         * qbasis (basisn, getcoeffxx)
                                     + (x - betax + 1) * (betax + 1)
                                           / ((t + 2) * (t + 1) * GG (0))
                                           * BB (x - betax + 1, y - betay)
                                           * qbasis (basisn, getcoeffx);
                              if (D == 2)
                                *newcoeff
                                    += (betay + 2) * (betay + 1)
                                           / ((t + 2) * (t + 1) * GG (0))
                                           * BB (x - betax, y - betay)
                                           * qbasis (basisn, getcoeffyy)
                                       + (y - betay + 1) * (betay + 1)
                                             / ((t + 2) * (t + 1) * GG (0))
                                             * BB (x - betax, y - betay + 1)
                                             * qbasis (basisn, getcoeffy);
                              if (betax + betay == x + y)
                                continue;
                              index[1] = betay;
                              index[0] = betax;
                              index[D] = t + 2;
                              int getcoeff = TrefftzWaveBasis<D>::IndexMap2 (
                                  index, ord);

                              *newcoeff -= GG (x - betax, y - betay)
                                           * qbasis (basisn, getcoeff)
                                           / GG (0);
                            }
                      }
                  }
              }
          }

        MatToCSR (qbasis, gtbstore[encode]);
      }

    if (gtbstore[encode].Size () == 0)
      {
        stringstream str;
        str << "failed to generate trefftz basis of order " << ord << endl;
        throw Exception (str.str ());
      }

    return gtbstore[encode];
  }

  template class QTrefftzWaveBasis<1>;
  template class QTrefftzWaveBasis<2>;
  template class QTrefftzWaveBasis<3>;

}
