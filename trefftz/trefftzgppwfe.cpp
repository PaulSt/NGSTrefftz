#include "trefftzgppwfe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
  template <int D>
  TrefftzGppwFE<D>::TrefftzGppwFE (int aord, float ac, Vec<D> aelcenter,
                                   double aelsize, ELEMENT_TYPE aeltype,
                                   int abasistype)
      : ScalarMappedElement<D + 1> (2 * aord + 1, aord), ord (aord), c (ac),
        npoly (BinCoeff (D + ord, ord)), elcenter (aelcenter),
        elsize (aelsize), eltype (aeltype), basistype (abasistype)
  {
    ;
  }

  template <>
  void TrefftzGppwFE<1>::CalcShape (const SIMD_BaseMappedIntegrationRule &smir,
                                    BareSliceMatrix<SIMD<double>> shape) const
  {
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
                = pow (GetDirection (i, d) * cpoint[0] - cpoint[1], i);
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
  }

  template <>
  void
  TrefftzGppwFE<1>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
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
                      * pow (GetDirection (i, dir) * cpoint[0] - cpoint[1],
                             (i - 1) * (i > 0))
                      * (d == 0 ? GetDirection (i, dir) : 1)
                      * (d == 1 ? (-c) : 1) * (2.0 / elsize);
          }
      }
  }

  template <>
  void
  TrefftzGppwFE<2>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
  }

  template <>
  void
  TrefftzGppwFE<3>::CalcDShape (const SIMD_BaseMappedIntegrationRule &smir,
                                BareSliceMatrix<SIMD<double>> dshape) const
  {
  }

  /////////////// non-simd

  template <>
  void TrefftzGppwFE<1>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<2>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
  }

  template <>
  void TrefftzGppwFE<3>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                    BareSliceVector<> shape) const
  {
  }

  template <>
  void TrefftzGppwFE<1>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
    cout << "dim not implemented" << endl;
  }

  template <>
  void TrefftzGppwFE<2>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
  }

  template <>
  void TrefftzGppwFE<3>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                     SliceMatrix<> dshape) const
  {
  }

  template class TrefftzGppwFE<1>;
  template class TrefftzGppwFE<2>;
  template class TrefftzGppwFE<3>;

  template <int D> const CSR *TrefftzGppwBasis<D>::TB (int ord)
  {
    const CSR *tb = &tbstore[ord];
    return tb;
  }

  template <int D>
  void TrefftzGppwBasis<D>::CreateTB (int ord, int gppword, Vector<> gamma,
                                      int basistype)
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
    const int nbasis = 2 * ord + 1;
    const int npoly = (gppword + 1) * (gppword + 2) / 2 - 1;
    Matrix<> trefftzbasis (nbasis, npoly);
    trefftzbasis = 0;
    Vec<D + 1, int> coeff = 0;
    int count = 0;
    int basisn = 0;
    for (int j = 0; j <= ord; ++j)
      for (int dir = 0; dir < TrefftzGppwFE<D>::NDirections (j); ++dir)
        for (int ell = (basisn++) - 1; ell < gppword; ell++)
          {
            // TB_inner(gamma, ord, trefftzbasis, coeff, b, D+1, ell,
            // basistype);
            Vec<D + 1, int> get_coeff;
            get_coeff[D] = 0;
            get_coeff[0] = ell + 2;
            trefftzbasis (basisn, IndexMap2 (get_coeff, gppword - 1)) = 0;
            get_coeff[D] = 1;
            get_coeff[0] = ell + 1;
            trefftzbasis (basisn, IndexMap2 (get_coeff, gppword - 1)) = 0;
            for (int t = 0; t < ell; t++)
              {
                int x = ell - t;
                get_coeff[D] = t + 2;
                get_coeff[0] = x;
                Vec<D + 1, int> get_coeff2;
                get_coeff2[D] = t;
                get_coeff2[0] = x + 2;
                trefftzbasis (basisn, IndexMap2 (get_coeff, gppword - 1))
                    = (x + 2) * (x + 1) / (t + 2) * (t + 2)
                          * trefftzbasis (basisn,
                                          IndexMap2 (get_coeff2, gppword - 1))
                      - (x <= basisn - 2) * BinCoeff (basisn, t);
              }
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
  void TrefftzGppwBasis<D>::TB_inner (Vector<> gamma, int ord,
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
