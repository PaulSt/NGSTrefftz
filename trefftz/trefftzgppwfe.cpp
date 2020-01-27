#include "trefftzgppwfe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"
#include "trefftzwavefe.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    TrefftzGppwFE<D> :: TrefftzGppwFE(const Array<double> &agamma, int aord, float ac, Vec<D+1> aelcenter, double aelsize, ELEMENT_TYPE aeltype, int abasistype)
    : ScalarMappedElement<D+1>(BinCoeff(D + aord, aord) + BinCoeff(D + aord-1, aord-1), aord),
    ord(aord),
    c(ac),
    npoly(BinCoeff(D+1 + ord, ord)),
    elcenter(aelcenter),
    elsize(aelsize),
    eltype(aeltype),
    basistype(abasistype),
    gamma(agamma)
    {;}


    template<>
    void TrefftzGppwFE<1> :: CalcShape (const SIMD_BaseMappedIntegrationRule & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    { 

        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize);
            Array<double> gam(gamma);
            gam[0] += elcenter[0];
            gam[1] *= (elsize/2.0);

            // calc 1 dimensional monomial basis
            STACK_ARRAY(SIMD<double>, mem, 2*(ord+1));
            Vec<2,SIMD<double>*> polxt;
            for(size_t d=0;d<2;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }
            // calc D+1 dimenional monomial basis
            Vector<SIMD<double>> pol(npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    pol[ii++] = polxt[0][i] * polxt[1][j];
            // TB*monomials for trefftz shape fcts
            const CSR* localmat = TrefftzGppwBasis<1>::getInstance().TB(ord,gam);
            for (int i=0; i<this->ndof; ++i)
            {
                shape(i,imip) = 0.0;
                for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                    shape(i,imip) += (*localmat)[2][j]*pol[(*localmat)[1][j]];
            }
        }
    }

    template<>
    void TrefftzGppwFE<2> :: CalcShape (const SIMD_BaseMappedIntegrationRule & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    { 
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<3,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize);
            Array<double> gam(gamma);
            gam[0] += elcenter[0];
            gam[1] *= (elsize/2.0);

            // calc 1 dimensional monomial basis
            STACK_ARRAY(SIMD<double>, mem, 3*(ord+1));
            Vec<3,SIMD<double>*> polxt;
            for(size_t d=0;d<3;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }
            // calc D+1 dimenional monomial basis
            Vector<SIMD<double>> pol(npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    for (size_t k = 0; k <= ord-i-j; k++)
                        pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];
            // TB*monomials for trefftz shape fcts
            const CSR* localmat = TrefftzGppwBasis<2>::getInstance().TB(ord,gam);
            for (int i=0; i<this->ndof; ++i)
            {
                shape(i,imip) = 0.0;
                for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                    shape(i,imip) += (*localmat)[2][j]*pol[(*localmat)[1][j]];
            }
        }
    }

    template<>
    void TrefftzGppwFE<3> :: CalcShape (const SIMD_BaseMappedIntegrationRule & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    { throw ExceptionNOSIMD("SIMD - CalcShape not overloaded"); }



    template<>
    void TrefftzGppwFE<1> :: CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
                                         BareSliceMatrix<SIMD<double>> dshape) const
    { 
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize);
            Array<double> gam(gamma);
            gam[0] += elcenter[0];
            gam[1] *= (elsize/2.0);

            // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
            STACK_ARRAY(SIMD<double>, mem, 2*(ord+1)+1); mem[0]=0;
            Vec<2,SIMD<double>*> polxt;
            for(size_t d=0;d<2;d++)
            {
                polxt[d] = &mem[d*(ord+1)+1];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            for(int d=0;d<2;d++)
            {
                Vector<SIMD<double>> pol(npoly);
                for (size_t i = 0, ii = 0; i <=ord; i++)
                    for (size_t j = 0; j <= ord-i; j++)
                        pol[ii++] = (d==0?i:(d==1?j:0))
                            * polxt[0][i-(d==0)] * polxt[1][j-(d==1)];

                const CSR* localmat = TrefftzGppwBasis<1>::getInstance().TB(ord,gam);
                for (int i=0; i<this->ndof; ++i)
                {
                    dshape(i*2+d,imip) = 0.0;
                    for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                        dshape(i*2+d,imip) += (*localmat)[2][j]*pol[(*localmat)[1][j]] * (d==1 ? c : 1) * (2.0/elsize);
                }
            }
        }
    }


    template<>
    void TrefftzGppwFE<2> :: CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
                                         BareSliceMatrix<SIMD<double>> dshape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<3,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize);
            Array<double> gam(gamma);
            gam[0] += elcenter[0];
            gam[1] *= (elsize/2.0);

            // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
            STACK_ARRAY(SIMD<double>, mem, 3*(ord+1)+1); mem[0]=0;
            Vec<3,SIMD<double>*> polxt;
            for(size_t d=0;d<3;d++)
            {
                polxt[d] = &mem[d*(ord+1)+1];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            for(int d=0;d<3;d++)
            {
                Vector<SIMD<double>> pol(npoly);
                for (size_t i = 0, ii = 0; i <=ord; i++)
                    for (size_t j = 0; j <= ord-i; j++)
                        for (size_t k = 0; k <= ord-i-j; k++)
                            pol[ii++] = (d==0?i:(d==1?j:(d==2?k:0)))
                                * polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * polxt[2][k-(d==2)];

                const CSR* localmat = TrefftzGppwBasis<2>::getInstance().TB(ord,gam);
                for (int i=0; i<this->ndof; ++i)
                {
                    dshape(i*3+d,imip) = 0.0;
                    for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                        dshape(i*3+d,imip) += (*localmat)[2][j]*pol[(*localmat)[1][j]] * (2.0/elsize);
                }
            }
        }
        //dshape *= (2.0/elsize); //inner derivative
    }
    template<>
    void TrefftzGppwFE<3> :: CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
                                         BareSliceMatrix<SIMD<double>> dshape) const
    { throw ExceptionNOSIMD("SIMD - CalcShape not overloaded"); }


    /////////////// non-simd


    template<>
    void TrefftzGppwFE<1> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
        Vec<2> cpoint = mip.GetPoint();
        cpoint -= elcenter;
        cpoint *= (2.0/elsize);
        Array<double> gam(gamma);
        gam[0] += elcenter[0];
        gam[1] *= (elsize/2.0);

        // calc 1 dimensional monomial basis
        STACK_ARRAY(double, mem, 2*(ord+1));
        double* polxt[2];
        for(size_t d=0;d<2;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }
        // calc D+1 dimenional monomial basis
        Vector<double> pol(npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
            for (size_t j = 0; j <= ord-i; j++)
                pol[ii++] = polxt[0][i] * polxt[1][j];
        // TB*monomials for trefftz shape fcts
        const CSR* localmat = TrefftzGppwBasis<1>::getInstance().TB(ord,gam);
        for (int i=0; i<this->ndof; ++i)
        {
            shape(i) = 0.0;
            for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                shape(i) += (*localmat)[2][j]*pol[(*localmat)[1][j]];
        }
    }

    template<>
    void TrefftzGppwFE<2> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
        Vec<3> cpoint = mip.GetPoint();
        cpoint -= elcenter;
        cpoint *= (2.0/elsize);
        Array<double> gam(gamma);
        gam[0] += elcenter[0]+elcenter[1];
        gam[1] *= (elsize/2.0);

        // calc 1 dimensional monomial basis
        STACK_ARRAY(double, mem, 3*(ord+1));
        double* polxt[3];
        for(size_t d=0;d<3;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }
        // calc D+1 dimenional monomial basis
        Vector<double> pol(npoly);
        for (size_t i = 0, ii = 0; i <= ord; i++)
            for (size_t j = 0; j <= ord-i; j++)
                for (size_t k = 0; k <= ord-i-j; k++)
                    pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];
        // TB*monomials for trefftz shape fcts
        const CSR* localmat = TrefftzGppwBasis<2>::getInstance().TB(ord,gam);
        for (int i=0; i<this->ndof; ++i)
        {
            shape(i) = 0.0;
            for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                shape(i) += (*localmat)[2][j]*pol[(*localmat)[1][j]];
        }
    }

    template<>
    void TrefftzGppwFE<3> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {cout << "dim not implemented" << endl;}



    template<>
    void TrefftzGppwFE<1> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         BareSliceMatrix<> dshape) const
    {
        Vec<2> cpoint = mip.GetPoint();
        cpoint -= elcenter;
        cpoint *= (2.0/elsize);
        Array<double> gam(gamma);
        gam[0] += elcenter[0];
        gam[1] *= (elsize/2.0);

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY(double, mem2, 2*(ord+1)+1); mem2[0]=0;
        int npoly = BinCoeff(1+1 + ord, ord);
        double* polxt2[2];
        for(size_t d=0;d<2;d++)
        {
            polxt2[d] = &mem2[d*(ord+1)+1];
            Monomial (ord, cpoint[d], polxt2[d]);
        }

        for(int d=0;d<2;d++)
        {
            Vector<double> pol(npoly);
            for (size_t i = 0, ii = 0; i <=ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    pol[ii++] = (d==0?i:(d==1?j:0))
                        * polxt2[0][i-(d==0)] * polxt2[1][j-(d==1)];

            const CSR* localmat = TrefftzGppwBasis<1>::getInstance().TB(ord,gam);
            for (int i=0; i<this->ndof; ++i)
            {
                dshape(i,d) = 0.0;
                for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                    dshape(i,d) += (*localmat)[2][j]*pol[(*localmat)[1][j]] * (2.0/elsize);
            }
        }
    }

    template<>
    void TrefftzGppwFE<2> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         BareSliceMatrix<> dshape) const
    {
        Vec<3> cpoint = mip.GetPoint();
        cpoint -= elcenter;
        cpoint *= (2.0/elsize);
        Array<double> gam(gamma);
        gam[0] += elcenter[0]+elcenter[1];
        gam[1] *= (elsize/2.0);

        // +1 size to avoid undefined behavior taking deriv, getting [-1] entry
        STACK_ARRAY(double, mem, 3*(ord+1)+1); mem[0]=0;
        double* polxt[3];
        for(size_t d=0;d<3;d++)
        {
            polxt[d] = &mem[d*(ord+1)+1];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        for(int d=0;d<3;d++)
        {
            Vector<double> pol(npoly);
            for (size_t i = 0, ii = 0; i <=ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    for (size_t k = 0; k <= ord-i-j; k++)
                        pol[ii++] = (d==0?i:(d==1?j:(d==2?k:0)))
                            * polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * polxt[2][k-(d==2)];

            const CSR* localmat = TrefftzGppwBasis<2>::getInstance().TB(ord,gam);
            for (int i=0; i<this->ndof; ++i)
            {
                dshape(i,d) = 0.0;
                for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                    dshape(i,d) += (*localmat)[2][j]*pol[(*localmat)[1][j]] * (2.0/elsize);
            }
        }
    }

    template<>
    void TrefftzGppwFE<3> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         BareSliceMatrix<> dshape) const
    {cout << "dim not implemented" << endl;}



    template class TrefftzGppwFE<1>;
    template class TrefftzGppwFE<2>;
    template class TrefftzGppwFE<3>;


    template<int D>
    const CSR* TrefftzGppwBasis<D> :: TB(int ord, const Array<double> &gamma, int basistype)
    {
        {
            lock_guard<mutex> lock(gentrefftzbasis);
            string encode = to_string(ord);
            for(auto g : gamma)
                encode += to_string(g);

            if ( gtbstore[encode][0].Size() == 0)
            {
                //cout << "creating gppw bstore for " << encode << endl;
                const int nbasis = (BinCoeff(D + ord, ord) + BinCoeff(D + ord-1, ord-1));
                const int npoly = BinCoeff(D+1 + ord, ord);
                Matrix<> gppwbasis(nbasis,npoly);
                gppwbasis = 0;

                Matrix<> trefftzbasis(nbasis,npoly);
                trefftzbasis = 0;
                Vec<D+1, int>  coeff = 0;
                int count = 0;
                for(int b=0;b<nbasis;b++)
                {
                    int tracker = 0;
                    TrefftzWaveBasis<D>::TB_inner(ord, trefftzbasis, coeff, b, D+1, tracker, basistype, 1/sqrt(gamma[0]));
                }

                for(int basisn=0;basisn<nbasis;basisn++)
                {
                    int j=0; // order of current basis fct
                    for (size_t i = 0; i <=ord; i++)
                        for (size_t k = 0; k <= (D==2)*(ord-i); k++)
                            for (size_t l = 0; l <= ord-i-k; l++)
                            {
                                Vec<D+1, int> index;
                                index[D] = l;
                                index[0] = i;
                                if(D==2) index[1] = k;
                                if (trefftzbasis( basisn, TrefftzWaveBasis<D>::IndexMap2(index, ord))!=0 && i+k+l>j)
                                    j=i+k+l;
                            }

                    for(int ell=-1;ell<ord-1;ell++)
                    {
                        for(int t=0;t<=ell;t++)
                        {
                            if(D==1)
                            {
                                int x = ell - t;
                                Vec<D+1, int> index;
                                index[D] = t+2;
                                index[0] = x;
                                double* newcoeff =& gppwbasis( basisn, TrefftzWaveBasis<D>::IndexMap2(index, ord));
                                index[D] = t;
                                index[0] = x+2;
                                int getcoeff = TrefftzWaveBasis<D>::IndexMap2(index, ord);

                                *newcoeff =
                                    (x+2)*(x+1)/((t+2)*(t+1)*gamma[0])
                                    * gppwbasis( basisn, getcoeff);
                                for(int betax=0;betax<x;betax++)
                                {
                                    index[D] = t+2;
                                    index[0] = betax;
                                    getcoeff = TrefftzWaveBasis<D>::IndexMap2(index, ord);

                                    *newcoeff
                                        -= gamma[x-betax]*gppwbasis( basisn, getcoeff) / gamma[0];
                                    if(t<=j-2)
                                        *newcoeff
                                            -= gamma[x-betax]*trefftzbasis( basisn, getcoeff) / gamma[0];
                                }
                            }
                            else if (D==2)
                            {
                                for(int x=0;x<=ell-t;x++)
                                {
                                    int y = ell-t-x;
                                    Vec<D+1, int> index;
                                    index[D] = t+2;
                                    index[1] = y;
                                    index[0] = x;
                                    double* newcoeff =& gppwbasis( basisn, TrefftzWaveBasis<D>::IndexMap2(index, ord));
                                    index[D] = t;
                                    index[1] = y;
                                    index[0] = x+2;
                                    int getcoeffx = TrefftzWaveBasis<D>::IndexMap2(index, ord);
                                    index[D] = t;
                                    index[1] = y+2;
                                    index[0] = x;
                                    int getcoeffy = TrefftzWaveBasis<D>::IndexMap2(index, ord);

                                    *newcoeff =
                                        (x+2)*(x+1)/((t+2)*(t+1)*gamma[0])
                                        * gppwbasis( basisn, getcoeffx)
                                        + (y+2)*(y+1)/((t+2)*(t+1)*gamma[0])
                                        * gppwbasis( basisn, getcoeffy);
                                    for(int betax=0;betax<=x;betax++)
                                        for(int betay=0;betay<=y-(betax==x);betay++)
                                        {
                                            index[D] = t+2;
                                            index[1] = betay;
                                            index[0] = betax;
                                            int getcoeff = TrefftzWaveBasis<D>::IndexMap2(index, ord);

                                            // TODO fix smart gamma, no only dummy
                                            double fakegamma = 0;
                                            if( (x-betax == 0 && y-betay == 0))
                                                fakegamma=gamma[0];
                                            else if((x-betax == 0 && y-betay == 1) || (x-betax == 1 && y-betay == 0))
                                                fakegamma=gamma[1];

                                            *newcoeff
                                                -= fakegamma*gppwbasis( basisn, getcoeff) / gamma[0];
                                            if(t<=j-2)
                                                *newcoeff
                                                    -= fakegamma*trefftzbasis( basisn, getcoeff) / gamma[0];
                                        }
                                }
                            }
                        }
                    }
                }

                for(int basisn=0;basisn<nbasis;basisn++)
                    for(int polyn=0;polyn<npoly;polyn++)
                        gppwbasis(basisn,polyn)+=trefftzbasis(basisn,polyn);

                MatToCSR(gppwbasis,gtbstore[encode]);
            }

            if ( gtbstore[encode].Size() == 0)
            {
                stringstream str;
                str << "failed to generate trefftz basis of order " << ord << endl;
                throw Exception (str.str());
            }

            const CSR* tb =& gtbstore[encode];
            return tb;
        }
    }

    template class TrefftzGppwBasis<1>;
    template class TrefftzGppwBasis<2>;
    template class TrefftzGppwBasis<3>;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//#ifdef NGS_PYTHON
//void ExportTrefftzElement(py::module m)
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
