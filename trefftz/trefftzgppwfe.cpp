#include "trefftzgppwfe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
  template <int D>
  TrefftzGppwFE<D>::TrefftzGppwFE (const Array<double> &agamma, int agppword,
                                   int aord, float ac, Vec<D> aelcenter,
                                   double aelsize, ELEMENT_TYPE aeltype,
                                   int abasistype)
      : ScalarMappedElement<D + 1> (2 * aord + 1, aord), ord (aord), c (ac),
        npoly (BinCoeff (D + 1 + ord, ord)), elcenter (aelcenter),
        elsize (aelsize), eltype (aeltype), basistype (abasistype),
        gppword (agppword), gamma (agamma)
  {
    ;
  }

  template <>
  void TrefftzGppwFE<1>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
    throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[1] *= c;
        // calc 1 dimensional monomial basis
        int basisn = 0;
        for (int i = 0; i <= ord; ++i)
          for (int d = 0; d < NDirections (i); ++d)
            shape (basisn++, imip)
                = pow (GetDirection (i, d) * cpoint[0] + cpoint[1], i);

        // calc 1 dimensional monomial basis
        STACK_ARRAY (SIMD<double>, mem, 2 * (gppword + 1));
        int npoly = BinCoeff (1 + 1 + gppword, gppword);
        Vec<2, SIMD<double> *> polxt;
        for (size_t d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (gppword + 1)];
            Monomial (gppword, cpoint[d], polxt[d]);
          }
        // calc D+1 dimenional monomial basis
        Vector<SIMD<double>> pol (npoly);
        for (size_t i = 0, ii = 0; i <= gppword; i++)
          for (size_t j = 0; j <= gppword - i; j++)
            pol[ii++] = polxt[0][i] * polxt[1][j];
        // TB*monomials for trefftz shape fcts
        const CSR *localmat
            = TrefftzGppwBasis<1>::getInstance ().TB (ord, gppword, gamma);
        for (int i = 0; i < this->ndof; ++i)
          {
            // shape(i,imip) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              shape (i, imip) += (*localmat)[2][j] * pol[(*localmat)[1][j]];
          }
      }
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
    for (int imip = 0; imip < smir.Size (); imip++)
      {
        Vec<2, SIMD<double>> cpoint = smir[imip].GetPoint ();
        cpoint -= elcenter;
        cpoint *= (2.0 / elsize);
        cpoint[1] *= c;
        // calc 1 dimensional monomial basis
        for (int d = 0; d < 2; d++)
          {
            int basisn = 0;
            for (int i = 0; i <= ord; ++i)
              for (int dir = 0; dir < NDirections (i); ++dir)
                dshape (2 * (basisn++) + d, imip)
                    = i
                      * pow (GetDirection (i, dir) * cpoint[0] + cpoint[1],
                             (i - 1) * (i > 0))
                      * (d == 0 ? GetDirection (i, dir) : 1)
                      * (d == 1 ? (c) : 1) * (2.0 / elsize);
          }

        //+1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY (SIMD<double>, mem, 2 * (gppword + 1) + 1);
        mem[0] = 0;
        int npoly = BinCoeff (1 + 1 + gppword, gppword);
        Vec<2, SIMD<double> *> polxt;
        for (size_t d = 0; d < 2; d++)
          {
            polxt[d] = &mem[d * (gppword + 1) + 1];
            Monomial (gppword, cpoint[d], polxt[d]);
          }

        for (int d = 0; d < 2; d++)
          {
            Vector<SIMD<double>> pol (npoly);
            for (size_t i = 0, ii = 0; i <= gppword; i++)
              for (size_t j = 0; j <= gppword - i; j++)
                pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                            * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

            const CSR *localmat
                = TrefftzGppwBasis<1>::getInstance ().TB (ord, gppword, gamma);
            for (int i = 0; i < this->ndof; ++i)
              {
                // dshape(i*2+d,imip) = 0.0;
                for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
                  dshape (i * 2 + d, imip)
                      += (*localmat)[2][j] * pol[(*localmat)[1][j]]
                         * (d == 1 ? c : 1) * (2.0 / elsize);
              }
          }
      }
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
    // Vec<2,double> cpoint = mip.GetPoint();
    // cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;
    //// calc 1 dimensional monomial basis
    // STACK_ARRAY(double, mem, 2*(ord+1)+1);
    // for(size_t d=0;d<2;d++)
    //{
    // double evalpoint = pow(-1,d)*cpoint[0]-cpoint[1];
    // Monomial (ord, evalpoint, &mem[d*(ord+1)]);
    //}
    // for (int i=0; i<this->ndof; ++i)
    //{
    // shape(i) = i<=ord ? mem[i] : mem[i+1];
    //}
    // shape(0)=1;
    // shape(1)=cpoint[0]-c*cpoint[1];
    // shape(2)=-cpoint[0]-c*cpoint[1];
    // int basisn=3;
    // for (int i=2; i<=ord;i++)
    // for(int d=0;d<NDirections(i);d++)
    // shape(basisn++)=shape(basisn-2)*((d%2?-1:1)*cpoint[0]-c*cpoint[1]);
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    // calc 1 dimensional monomial basis
    for (int i = 0, basisn = 0; i <= ord; ++i)
      for (int d = 0; d < NDirections (i); ++d)
        shape (basisn++)
            = pow (GetDirection (i, d) * cpoint[0] + c * cpoint[1], i);

    // calc 1 dimensional monomial basis
    STACK_ARRAY (double, mem, 2 * (gppword + 1));
    int npoly = BinCoeff (1 + 1 + gppword, gppword);
    double *polxt[2];
    for (size_t d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (gppword + 1)];
        Monomial (gppword, cpoint[d], polxt[d]);
      }
    // calc D+1 dimenional monomial basis
    Vector<double> pol (npoly);
    for (size_t i = 0, ii = 0; i <= gppword; i++)
      for (size_t j = 0; j <= gppword - i; j++)
        pol[ii++] = polxt[0][i] * polxt[1][j];
    // TB*monomials for trefftz shape fcts
    const CSR *localmat
        = TrefftzGppwBasis<1>::getInstance ().TB (ord, gppword, gamma);
    for (int i = 0; i < this->ndof; ++i)
      {
        // shape(i) = 0.0;
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
                                     SliceMatrix<> dshape) const
  {
    Vec<2> cpoint = mip.GetPoint ();
    cpoint -= elcenter;
    cpoint *= (2.0 / elsize);
    // calc 1 dimensional monomial basis
    for (int d = 0; d < 2; d++)
      {
        for (int i = 0, basisn = 0; i <= ord; ++i)
          for (int dir = 0; dir < NDirections (i); ++dir)
            dshape (basisn++, d)
                = i
                  * pow (GetDirection (i, dir) * cpoint[0] + c * cpoint[1],
                         (i - 1) * (i > 0))
                  * (d == 0 ? GetDirection (i, dir) : 1) * (d == 1 ? (c) : 1)
                  * (2.0 / elsize);
      }

    // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
    STACK_ARRAY (double, mem, 2 * (gppword + 1) + 1);
    mem[0] = 0;
    int npoly = BinCoeff (1 + 1 + gppword, gppword);
    double *polxt[2];
    for (size_t d = 0; d < 2; d++)
      {
        polxt[d] = &mem[d * (gppword + 1) + 1];
        Monomial (gppword, cpoint[d], polxt[d]);
      }

    for (int d = 0; d < 2; d++)
      {
        Vector<double> pol (npoly);
        for (size_t i = 0, ii = 0; i <= gppword; i++)
          for (size_t j = 0; j <= gppword - i; j++)
            pol[ii++] = (d == 0 ? i : (d == 1 ? j : 0))
                        * polxt[0][i - (d == 0)] * polxt[1][j - (d == 1)];

        const CSR *localmat
            = TrefftzGppwBasis<1>::getInstance ().TB (ord, gppword, gamma);
        for (int i = 0; i < this->ndof; ++i)
          {
            // dshape(i,d) = 0.0;
            for (int j = (*localmat)[0][i]; j < (*localmat)[0][i + 1]; ++j)
              dshape (i, d) += (*localmat)[2][j] * pol[(*localmat)[1][j]]
                               * (2.0 / elsize);
          }
      }
  }

  template <>
  void TrefftzGppwFE<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<3>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template class TrefftzGppwFE<1>;
  template class TrefftzGppwFE<2>;
  template class TrefftzGppwFE<3>;

  template <int D>
  const CSR *
  TrefftzGppwBasis<D>::TB (int ord, int gppword, const Array<double> &gamma,
                           int basistype)
  {
    {
      lock_guard<mutex> lock (gentrefftzbasis);
      string encode = to_string (ord);
      for (auto g : gamma)
        encode += to_string (g);

      if (gtbstore[encode][0].Size () == 0)
        {
          cout << "creating gppw bstore for " << encode << endl;
          const int nbasis = 2 * ord + 1;
          const int npoly = BinCoeff (D + 1 + gppword, gppword);
          Matrix<> trefftzbasis (nbasis, npoly);
          trefftzbasis = 0;
          Vec<D + 1, int> coeff = 0;
          int count = 0;
          int basisn = 0;
          for (int j = 0; j <= ord; j++)
            {
              for (int dir = 0; dir < TrefftzGppwFE<D>::NDirections (j); dir++)
                {
                  for (int ell = j - 1; ell <= gppword - 2; ell++)
                    {
                      Vec<D + 1, int> get_coeff;
                      get_coeff[D] = 0;
                      get_coeff[0] = ell + 2;
                      trefftzbasis (basisn, IndexMap2 (get_coeff, gppword))
                          = 0;
                      get_coeff[D] = 1;
                      get_coeff[0] = ell + 1;
                      trefftzbasis (basisn, IndexMap2 (get_coeff, gppword))
                          = 0;
                      for (int t = 0; t <= ell; t++)
                        {
                          int x = ell - t;
                          get_coeff[D] = t + 2;
                          get_coeff[0] = x;

                          // cout << "setting " << get_coeff << " mapped to "
                          // << IndexMap2(get_coeff, gppword) << endl;
                          Vec<D + 1, int> get_coeff2;
                          get_coeff2[D] = t;
                          get_coeff2[0] = x + 2;
                          trefftzbasis (basisn, IndexMap2 (get_coeff, gppword))
                              = (x + 2) * (x + 1)
                                    / ((t + 2) * (t + 1) * gamma[0])
                                    * trefftzbasis (
                                        basisn,
                                        IndexMap2 (get_coeff2, gppword))
                                - (t <= j - 2) * BinCoeff (j, t + 2)
                                      * gamma[x + t + 2 - j]
                                      * pow (gamma[0], 0.5 * (j - t - 4))
                                      * pow (TrefftzGppwFE<D>::GetDirection (
                                                 j, dir),
                                             x);
                          // if(t<=j-2) cout << "ell " << ell << " ord " << j
                          // << " at " << t << endl;
                          for (int betax = max (0, j - t - 1); betax < x;
                               betax++)
                            {
                              get_coeff2[D] = t + 2;
                              get_coeff2[0] = betax;
                              trefftzbasis (basisn,
                                            IndexMap2 (get_coeff, gppword))
                                  -= gamma[x - betax]
                                     * trefftzbasis (
                                         basisn,
                                         IndexMap2 (get_coeff2, gppword))
                                     / gamma[0];
                            }
                        }
                    }
                  basisn++;
                }
            }
          cout << trefftzbasis << endl;
          cout << "size " << trefftzbasis.Height () << " x "
               << trefftzbasis.Width () << endl;

          MatToCSR (trefftzbasis, gtbstore[encode]);
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

  template <int D>
  void TrefftzGppwBasis<D>::TB_inner (const Array<double> &gamma, int ord,
                                      Matrix<> &trefftzbasis,
                                      Vec<D + 1, int> coeffnum, int basis,
                                      int dim, int &tracker, int basistype)
  {
    // if (dim>0)
    //{
    // while(coeffnum(dim-1)<=ord)
    //{
    // TB_inner(gamma, ord,trefftzbasis,coeffnum,basis, dim-1, tracker,
    // basistype); coeffnum(dim-1)++;
    // }
    // }
    // else
    //{
    // int sum=0;
    // for(int i=0;i<D+1;i++)
    // sum += coeffnum(i);
    // if(sum<=ord)
    //{
    // int indexmap = IndexMap2(coeffnum, ord);
    // int k = coeffnum(D);
    // cout << coeffnum << endl;
    ////trefftzbasis(basis, indexmap) *= 1.0/(k * (k-1));
    //}
    //}
  }

  template <int D>
  int TrefftzGppwBasis<D>::IndexMap2 (Vec<D + 1, int> index, int ord)
  {
    int sum = 0;
    int temp_size = 0;
    for (int d = 0; d < D + 1; d++)
      {
        for (int p = 0; p < index (d); p++)
          {
            sum += BinCoeff (D - d + ord - p - temp_size, ord - p - temp_size);
          }
        temp_size += index (d);
      }
    return sum;
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
