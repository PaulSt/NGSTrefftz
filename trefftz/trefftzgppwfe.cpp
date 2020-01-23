#include "trefftzgppwfe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"
#include "trefftzwavefe.hpp"

#include <ctime>

namespace ngfem
{
  template <int D>
  TrefftzGppwFE<D>::TrefftzGppwFE (const Array<double> &agamma, int aord,
                                   float ac, Vec<D + 1> aelcenter,
                                   double aelsize, ELEMENT_TYPE aeltype,
                                   int abasistype)
      : ScalarMappedElement<D + 1> (
          BinCoeff (D + aord, aord) + BinCoeff (D + aord - 1, aord - 1), aord),
        ord (aord), c (ac), npoly (BinCoeff (D + 1 + ord, ord)),
        elcenter (aelcenter), elsize (aelsize), eltype (aeltype),
        basistype (abasistype), gamma (agamma)
  {
    ;
  }

  template <>
  void TrefftzGppwFE<1>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
  }

  template <>
  void TrefftzGppwFE<2>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<3>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void
  TrefftzGppwFE<1>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
    throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
  }

  template <>
  void
  TrefftzGppwFE<2>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void
  TrefftzGppwFE<3>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  /////////////// non-simd

  template <>
  void TrefftzGppwFE<1>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    Array<double> gam (gamma);
    gam[0] += elcenter[0];
    gam[1] *= (elsize / 2.0);

    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 2 * (ord + 1));
    int npoly = BinCoeff (1 + 1 + ord, ord);
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
    const CSR *localmat = TrefftzGppwBasis<1>::getInstance ().TB (ord, gam);
    for (int i = 0; i < this->ndof; ++i)
      {
        shape (i) = 0.0;
        for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
          shape (i) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
      }
  }

  template <>
  void TrefftzGppwFE<2>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<3>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<1>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceMatrix<> dshape) const
  {
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    Array<double> gam (gamma);
    gam[0] += elcenter[0];
    gam[1] *= (elsize / 2.0);

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem2, 2 * (ord + 1) + 1);
    mem2[0] = 0;
    int npoly = BinCoeff (1 + 1 + ord, ord);
    double *polxt2[2];
    for (size_t d = 0; d < 2; d++)
      {
        polxt2[d] = &mem2[d * (ord + 1) + 1];
        Monomial (ord, cpoint[d], polxt2[d]);
      }

    for (int d = 0; d < 2; d++)
      {
        Vector<double> pol (npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
          for (size_t j = 0; j <= ord - i; j++)
            pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                        * polxt2[0][i - (d == 0)] * polxt2[1][j - (d == 1)];

        const CSR *localmat
            = TrefftzGppwBasis<1>::getInstance ().TB (ord, gam);
        for (int i = 0; i < this->ndof; ++i)
          {
            dshape (i, d) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              dshape (i, d) += (*localmat)[2][j] * pol[(*localmat)[1][j]]
                               * (2.0 / elsize);
          }
      }
  }

  template <>
  void TrefftzGppwFE<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceMatrix<> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<3>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     BareSliceMatrix<> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template class TrefftzGppwFE<1>;
  template class TrefftzGppwFE<2>;
  template class TrefftzGppwFE<3>;

  template <int D>
  const CSR *
  TrefftzGppwBasis<D>::TB (int ord, const Array<double> &gamma, int basistype)
  {
    {
      lock_guard<mutex> lock (gentrefftzbasis);
      string encode = to_string (ord);
      for (auto g : gamma)
        encode += to_string (g);

      if (gtbstore[encode][0].Size () == 0)
        {
          // cout << "creating gppw bstore for " << encode << endl;
          const int nbasis
              = (BinCoeff (D + ord, ord) + BinCoeff (D + ord - 1, ord - 1));
          const int npoly = BinCoeff (D + 1 + ord, ord);
          Matrix<> gppwbasis (nbasis, npoly);
          gppwbasis = 0;

          Matrix<> trefftzbasis (nbasis, npoly);
          trefftzbasis = 0;
          Vec<D + 1, int> coeff = 0;
          int count = 0;
          for (int b = 0; b < nbasis; b++)
            {
              int tracker = 0;
              TrefftzWaveBasis<D>::TB_inner (ord, trefftzbasis, coeff, b,
                                             D + 1, tracker, basistype,
                                             1 / sqrt (gamma[0]));
            }

          for (int basisn = 0; basisn < nbasis; basisn++)
            {
              Vec<D + 1, int> get_coeff;
              int j = 0; // order of current basis fct
              for (size_t i = 0; i <= ord; i++)
                {
                  for (size_t k = 0; k <= ord - i; k++)
                    {
                      get_coeff[D] = k;
                      get_coeff[0] = i;
                      if (trefftzbasis (
                              basisn,
                              TrefftzWaveBasis<D>::IndexMap2 (get_coeff, ord))
                              != 0
                          && i + k > j)
                        j = i + k;
                    }
                }

              for (int ell = -1; ell < ord - 1; ell++)
                {
                  Vec<D + 1, int> get_coeff;
                  get_coeff[D] = 0;
                  get_coeff[0] = ell + 2;
                  gppwbasis (basisn,
                             TrefftzWaveBasis<D>::IndexMap2 (get_coeff, ord))
                      = 0;
                  get_coeff[D] = 1;
                  get_coeff[0] = ell + 1;
                  gppwbasis (basisn,
                             TrefftzWaveBasis<D>::IndexMap2 (get_coeff, ord))
                      = 0;
                  for (int t = 0; t <= ell; t++)
                    {
                      int x = ell - t;
                      get_coeff[D] = t + 2;
                      get_coeff[0] = x;
                      Vec<D + 1, int> get_coeff2;
                      get_coeff2[D] = t;
                      get_coeff2[0] = x + 2;

                      gppwbasis (basisn, TrefftzWaveBasis<D>::IndexMap2 (
                                             get_coeff, ord))
                          = (x + 2) * (x + 1) / ((t + 2) * (t + 1) * gamma[0])
                            * gppwbasis (basisn,
                                         TrefftzWaveBasis<D>::IndexMap2 (
                                             get_coeff2, ord));
                      for (int betax = 0; betax < x; betax++)
                        {
                          get_coeff2[D] = t + 2;
                          get_coeff2[0] = betax;

                          gppwbasis (basisn, TrefftzWaveBasis<D>::IndexMap2 (
                                                 get_coeff, ord))
                              -= gamma[x - betax]
                                 * gppwbasis (basisn,
                                              TrefftzWaveBasis<D>::IndexMap2 (
                                                  get_coeff2, ord))
                                 / gamma[0];
                          if (t <= j - 2)
                            gppwbasis (basisn, TrefftzWaveBasis<D>::IndexMap2 (
                                                   get_coeff, ord))
                                -= gamma[x - betax]
                                   * trefftzbasis (
                                       basisn, TrefftzWaveBasis<D>::IndexMap2 (
                                                   get_coeff2, ord))
                                   / gamma[0];
                        }
                    }
                }
            }

          for (int basisn = 0; basisn < nbasis; basisn++)
            for (int polyn = 0; polyn < npoly; polyn++)
              gppwbasis (basisn, polyn) += trefftzbasis (basisn, polyn);

          MatToCSR (gppwbasis, gtbstore[encode]);
        }

      if (gtbstore[encode].Size () == 0)
        {
          stringstream str;
          str << "failed to generate trefftz basis of order " << ord << endl;
          throw Exception (str.str ());
        }

      const CSR *tb = &gtbstore[encode];
      return tb;
    }
  }

  template class TrefftzGppwBasis<1>;
  template class TrefftzGppwBasis<2>;
  template class TrefftzGppwBasis<3>;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#ifdef NGS_PYTHON
// void ExportTrefftzElement(py::module m)
//{
//// py::class_<TrefftzGppwFE<3>, shared_ptr<TrefftzGppwFE<3>>, FiniteElement>
//// 	(m, "TrefftzWaveFE3", "Trefftz space for wave eq")
//// 	.def(py::init<>())
//// 	;
//// py::class_<TrefftzGppwFE<2>, shared_ptr<TrefftzGppwFE<2>>, FiniteElement>
//// 	(m, "TrefftzWaveFE2", "Trefftz space for wave eq")
//// 	.def(py::init<>())
//// 	;
//// py::class_<TrefftzGppwFE<1>, shared_ptr<TrefftzGppwFE<1>>, FiniteElement>
//// 	(m, "TrefftzWaveFE1", "Trefftz space for wave eq")
//// 	.def(py::init<>())
//// 	;
//}
//#endif // NGS_PYTHON
