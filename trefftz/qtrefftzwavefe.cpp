#include "qtrefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"
#include "trefftzwavefe.hpp"

#include <ctime>

namespace ngfem
{
  constexpr int factorial (int n) { return n > 1 ? n * factorial (n - 1) : 1; }

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
    string encode = to_string (ord) + to_string (elsize);
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

}
