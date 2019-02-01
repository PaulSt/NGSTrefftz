#include "trefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    TrefftzWaveFE<D> :: TrefftzWaveFE(int aord, float ac, Vec<D> aelcenter, double aelsize, ELEMENT_TYPE aeltype)
    : ScalarMappedElement<D>(BinCoeff(D-1 + aord, aord) + BinCoeff(D-1 + aord-1, aord-1), aord),
    ord(aord),
    c(ac),
    nbasis(BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1)),
    npoly(BinCoeff(D + ord, ord)),
    elcenter(aelcenter),
    elsize(aelsize),
    eltype(aeltype)
    {;}

    template<>
    void TrefftzWaveFE<1> :: CalcShape (const SIMD_MappedIntegrationRule<0,1> & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {cout << "dim not implemented" << endl;}

    template<>
    void TrefftzWaveFE<2> :: CalcShape (const SIMD_MappedIntegrationRule<1,2> & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;

            STACK_ARRAY(SIMD<double>, mem, 2*(ord+1));
            Vec<2,SIMD<double>*> polxt;
            for(size_t d=0;d<2;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            Vector<SIMD<double>> tempshape(nbasis);
            Vector<SIMD<double>> pol(npoly);

            for (size_t i = 0, ii = 0; i <= ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    pol[ii++] = polxt[0][i] * polxt[1][j];

            //Matrix<> localmat = *(TrefftzWaveBasis<2>::getInstance().TB(ord));
            //tempshape = localmat * pol;
            for(int b=0;b<nbasis;b++) shape.Col(imip)(b) = tempshape(b);
        }
    }

    template<>
    void TrefftzWaveFE<3> :: CalcShape (const SIMD_MappedIntegrationRule<2,3> & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<3,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[2] *= c;

            STACK_ARRAY(SIMD<double>, mem, 3*(ord+1));
            Vec<3,SIMD<double>*> polxt;
            for(size_t d=0;d<3;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            //Vector<SIMD<double>> tempshape(nbasis);
            Vector<SIMD<double>> pol(npoly);

            for (size_t i = 0, ii = 0; i <= ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    for (size_t k = 0; k <= ord-i-j; k++)
                        pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];

            //tempshape = 0.0;
            //for (int i=0; i<(*localmat)[0].Size(); i++)
                //tempshape[(*localmat)[0][i]] += (*localmat)[2][i]*pol[(*localmat)[1][i]];
            //for(int b=0;b<nbasis;b++) shape.Col(imip)(b) = tempshape(b);
            const CSR* localmat = TrefftzWaveBasis<3>::getInstance().TB(ord);

            for (int i=0; i<nbasis; ++i)
            {
                shape(i,imip) = 0.0;
                for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                    shape(i,imip) += (*localmat)[2][j]*pol[(*localmat)[1][j]];
            }
        }
    }

    template<>
    void TrefftzWaveFE<4> :: CalcShape (const SIMD_MappedIntegrationRule<3,4> & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<4,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[3] *= c;

            STACK_ARRAY(SIMD<double>, mem, 4*(ord+1));
            Vec<4,SIMD<double>*> polxt;
            for(size_t d=0;d<4;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            Vector<SIMD<double>> pol(npoly);
            for (size_t i = 0, ii = 0; i <= ord; i++)
                for (size_t j = 0; j <= ord-i; j++)
                    for (size_t k = 0; k <= ord-i-j; k++)
                        for (size_t l = 0; l <= ord-i-j-k; l++)
                            pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k] * polxt[3][l];

            const CSR* localmat = TrefftzWaveBasis<4>::getInstance().TB(ord);
            //for (int i=0; i<(*localmat)[0].Size(); i++)
            //shape(imip, (*localmat)[0][i]) += (*localmat)[2][i]*pol[(*localmat)[1][i]];

            for (int i=0; i<(*localmat)[0].Size()-1; ++i)
            {
                shape(i,imip) = 0.0;
                for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                    shape(i,imip) += (*localmat)[2][j]*pol[(*localmat)[1][j]];
            }
        }
    }

    template<>
    void TrefftzWaveFE<1> :: CalcDShape (const SIMD_MappedIntegrationRule<0,1> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {}

    template<>
    void TrefftzWaveFE<2> :: CalcDShape (const SIMD_MappedIntegrationRule<1,2> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;

            STACK_ARRAY(SIMD<double>, mem, 2*(ord+1));
            Vec<2,SIMD<double>*> polxt;
            for(size_t d=0;d<2;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            //Matrix<> localmat = *(TrefftzWaveBasis<2>::getInstance().TB(ord));
            //Vector<SIMD<double>> tempdshape(nbasis);
            //for(int d=0;d<2;d++)
            //{
            //tempdshape=0;
            //for (size_t i = 0, ii = 0; i <=ord; i++)
            //for (size_t j = 0; j <= ord-i; j++)
            //{
            //ii++;
            //if((d==0&&i==0)||(d==1&&j==0)) continue;
            //SIMD<double> pol = polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * (d==0?i:j);
            //tempdshape += pol * localmat.Col(ii-1);
            //}
            //for(int n=0;n<nbasis;n++)
            //dshape(n*2+d,imip) = tempdshape(n) * (d==1 ? c : 1);
            //}
        }
        dshape *= (2.0/elsize); //inner derivative
    }

    template<>
    void TrefftzWaveFE<3> :: CalcDShape (const SIMD_MappedIntegrationRule<2,3> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<3,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[2] *= c;

            STACK_ARRAY(SIMD<double>, mem, 3*(ord+1));
            Vec<3,SIMD<double>*> polxt;
            for(size_t d=0;d<3;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            //CSR localmat = *(TrefftzWaveBasis<3>::getInstance().TB(ord));

            for(int d=0;d<3;d++)
            {
                Vector<SIMD<double>> pol(npoly);
                for (size_t i = 0, ii = 0; i <=ord; i++)
                    for (size_t j = 0; j <= ord-i; j++)
                        for (size_t k = 0; k <= ord-i-j; k++)
                        {
                            if((d==0&&i==0)||(d==1&&j==0)||(d==2&&k==0))
                                pol[ii++] = 0;
                            else
                                pol[ii++] = polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * polxt[2][k-(d==2)] * (d==0?i:(d==1?j:k));
                        }
                Vector<SIMD<double>> tempdshape(nbasis);
                //tempdshape = 0.0;
                //for (int i=0; i<(*localmat)[0].Size(); i++)
                    //tempdshape[(*localmat)[0][i]] += (*localmat)[2][i]*pol[(*localmat)[1][i]];
                const CSR* localmat = TrefftzWaveBasis<3>::getInstance().TB(ord);
                for (int i=0; i<nbasis; ++i)
                {
                    tempdshape[i] = 0.0;
                    for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                        tempdshape[i] += (*localmat)[2][j]*pol[(*localmat)[1][j]];
                }
                for(int n=0;n<nbasis;n++)
                    dshape(n*3+d,imip) = tempdshape(n) * (d==2 ? c : 1);
            }
        }
        dshape *= (2.0/elsize); //inner derivative
    }

    template<>
    void TrefftzWaveFE<4> :: CalcDShape (const SIMD_MappedIntegrationRule<3,4> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<4,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[3] *= c;

            STACK_ARRAY(SIMD<double>, mem, 4*(ord+1));
            Vec<4,SIMD<double>*> polxt;
            for(size_t d=0;d<4;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            for(int d=0;d<4;d++)
            {
                Vector<SIMD<double>> pol(npoly);
                for (size_t i = 0, ii = 0; i <=ord; i++)
                    for (size_t j = 0; j <= ord-i; j++)
                        for (size_t k = 0; k <= ord-i-j; k++)
                            for (size_t l = 0; l <= ord-i-j-k; l++)
                            {
                                if((d==0&&i==0)||(d==1&&j==0)||(d==2&&k==0)||(d==3&&l==0))
                                    pol[ii++] = 0;
                                else
                                    pol[ii++] = polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * polxt[2][k-(d==2)] * polxt[3][l-(d==3)]
                                        * (d==0?i:(d==1?j:(d==2?k:l)));
                            }

                const CSR* localmat = TrefftzWaveBasis<4>::getInstance().TB(ord);
                Vector<SIMD<double>> tempdshape(nbasis);
                //for (int i=0; i<(*localmat)[0].Size(); i++)
                //tempdshape[(*localmat)[0][i]] += (*localmat)[2][i]*pol[(*localmat)[1][i]];
                for (int i=0; i<(*localmat)[0].Size()-1; ++i)
                {
                    tempdshape[i] = 0.0;
                    for (int j=(*localmat)[0][i]; j<(*localmat)[0][i+1]; ++j)
                        tempdshape[i] += (*localmat)[2][j]*pol[(*localmat)[1][j]];
                }
                for(int n=0;n<nbasis;n++)
                    dshape(n*4+d,imip) = tempdshape(n) * (d==3 ? c : 1);
            }
        }
        dshape *= (2.0/elsize); //inner derivative
    }

    /////////////// non-simd


    template<>
    void TrefftzWaveFE<1> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {cout << "dim not implemented" << endl;}

    template<>
    void TrefftzWaveFE<2> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        Vec<2> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;

        STACK_ARRAY(double, mem, 2*(ord+1));
        Vec<2,double*> polxt;
        for(size_t d=0;d<2;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        Vector<> tempshape(nbasis);
        Vector<> pol(npoly);

        for (size_t i = 0, ii = 0; i <= ord; i++)
            for (size_t j = 0; j <= ord-i; j++)
                pol[ii++] = polxt[0][i] * polxt[1][j];

        //Matrix<> localmat = *(TrefftzWaveBasis<2>::getInstance().TB(ord));
        //tempshape = localmat * pol;
        for(int b=0;b<nbasis;b++) shape(b) = tempshape(b);
    }

    template<>
    void TrefftzWaveFE<3> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        Vec<3> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[2] *= c;

        STACK_ARRAY(double, mem, 3*(ord+1));
        Vec<3,double*> polxt;
        for(size_t d=0;d<3;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        Vector<> tempshape(nbasis);
        Vector<> pol(npoly);

        for (size_t i = 0, ii = 0; i <= ord; i++)
            for (size_t j = 0; j <= ord-i; j++)
                for (size_t k = 0; k <= ord-i-j; k++)
                    pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k];

        //Matrix<> localmat = *(TrefftzWaveBasis<3>::getInstance().TB(ord));
        //tempshape = localmat * pol;
        for(int b=0;b<nbasis;b++) shape(b) = tempshape(b);
    }

    template<>
    void TrefftzWaveFE<4> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        Vec<4> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[3] *= c;

        STACK_ARRAY(double, mem, 4*(ord+1));
        Vec<4,double*> polxt;
        for(size_t d=0;d<4;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        Vector<> tempshape(nbasis);
        Vector<> pol(npoly);
        pol= 1;

        for (size_t i = 0, ii = 0; i <= ord; i++)
            for (size_t j = 0; j <= ord-i; j++)
                for (size_t k = 0; k <= ord-i-j; k++)
                    for (size_t l = 0; l <= ord-i-j-k; l++)
                        pol[ii++] = polxt[0][i] * polxt[1][j] * polxt[2][k] * polxt[3][l];

        //Matrix<> localmat = *(TrefftzWaveBasis<4>::getInstance().TB(ord));
        //tempshape = localmat * pol;
        for(int b=0;b<nbasis;b++) shape(b) = tempshape(b);
    }



    template<>
    void TrefftzWaveFE<1> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {}

    template<>
    void TrefftzWaveFE<2> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {
        Vec<2> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;

        STACK_ARRAY(double, mem, 2*(ord+1));
        Vec<2,double*> polxt;
        for(size_t d=0;d<2;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        //Matrix<> localmat = *(TrefftzWaveBasis<2>::getInstance().TB(ord));
        //Vector<> tempdshape(nbasis);
        //for(int d=0;d<2;d++)
        //{
        //tempdshape=0;
        //for (size_t i = 0, ii = 0; i <=ord; i++)
        //for (size_t j = 0; j <= ord-i; j++)
        //{
        //ii++;
        //if((d==0&&i==0)||(d==1&&j==0)) continue;
        //double pol = polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * (d==0?i:j);
        //tempdshape += pol * localmat.Col(ii-1);
        //}
        //dshape.Col(d) = tempdshape;
        //}
        dshape.Col(1) *= c; //inner derivative
        dshape *= (2.0/elsize); //inner derivative
    }

    template<>
    void TrefftzWaveFE<3> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {
        Vec<3> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[2] *= c;

        STACK_ARRAY(double, mem, 3*(ord+1));
        Vec<3,double*> polxt;
        for(size_t d=0;d<3;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        //Matrix<> localmat = *(TrefftzWaveBasis<3>::getInstance().TB(ord));
        //Vector<> tempdshape(nbasis);
        //for(int d=0;d<3;d++)
        //{
        //tempdshape=0;
        //for (size_t i = 0, ii = 0; i <=ord; i++)
        //for (size_t j = 0; j <= ord-i; j++)
        //for (size_t k = 0; k <= ord-i-j; k++)
        //{
        //ii++;
        //if((d==0&&i==0)||(d==1&&j==0)||(d==2&&k==0)) continue;
        //double pol = polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * polxt[2][k-(d==2)] * (d==0?i:(d==1?j:k));
        //tempdshape += pol * localmat.Col(ii-1);
        //}
        //dshape.Col(d) = tempdshape;
        //}
        dshape.Col(2) *= c; //inner derivative
        dshape *= (2.0/elsize); //inner derivative
    }

    template<>
    void TrefftzWaveFE<4> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {
        Vec<4> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[3] *= c;

        STACK_ARRAY(double, mem, 4*(ord+1));
        Vec<4,double*> polxt;
        for(size_t d=0;d<4;d++)
        {
            polxt[d] = &mem[d*(ord+1)];
            Monomial (ord, cpoint[d], polxt[d]);
        }

        //Matrix<> localmat = *(TrefftzWaveBasis<4>::getInstance().TB(ord));
        //Vector<> tempdshape(nbasis);
        //for(int d=0;d<4;d++)
        //{
        //tempdshape=0;
        //for (size_t i = 0, ii = 0; i <=ord; i++)
        //for (size_t j = 0; j <= ord-i; j++)
        //for (size_t k = 0; k <= ord-i-j; k++)
        //for (size_t l = 0; l <= ord-i-j-k; l++)
        //{
        //ii++;
        //if((d==0&&i==0)||(d==1&&j==0)||(d==2&&k==0)||(d==3&&l==0)) continue;
        //double pol = polxt[0][i-(d==0)] * polxt[1][j-(d==1)] * polxt[2][k-(d==2)] * polxt[3][l-(d==3)]
        //* (d==0?i:(d==1?j:(d==2?k:l)));
        //tempdshape += pol * localmat.Col(ii-1);
        //}
        //dshape.Col(d) = tempdshape;
        //}
        dshape.Col(3) *= c; //inner derivative
        dshape *= (2.0/elsize); //inner derivative
    }

    template class TrefftzWaveFE<1>;
    template class TrefftzWaveFE<2>;
    template class TrefftzWaveFE<3>;
    template class TrefftzWaveFE<4>;


    template<int D>
    const CSR* TrefftzWaveBasis<D> :: TB(int ord)
    {
        {
            lock_guard<mutex> lock(gentrefftzbasis);
            if (tbstore.Size() <= ord)
            {
                int oldsize = tbstore.Size();
                tbstore.SetSize (ord+1);
                for (int i = oldsize; i <= ord; i++)
                    tbstore[i] = CSR();
            }

            if ( tbstore[ord][0].Size() == 0)
            {
                cout << "basis for ord: " << ord << endl;
                const int nbasis = (BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1));
                const int npoly = (BinCoeff(D + ord, ord));
                Matrix<> trefftzbasis(nbasis,npoly);
                trefftzbasis = 0;
                Vec<D, int>  coeff = 0;
                int count = 0;
                for(int b=0;b<nbasis;b++)
                {
                    int tracker = 0;
                    TB_inner(ord, trefftzbasis, coeff, b, D, tracker);
                }
                MatToCOO(trefftzbasis,tbstore[ord]);
            }

            if ( tbstore[ord][0].Size() == 0)
            {
                stringstream str;
                str << "failed to generate trefftz basis of order " << ord << endl;
                throw Exception (str.str());
            }
        }
        const CSR* tb =& tbstore[ord];
        return tb;
    }


    template<int D>
    void TrefftzWaveBasis<D> :: TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D, int> coeffnum, int basis, int dim, int &tracker)
    {
        if (dim>0)
        {
            while(coeffnum(dim-1)<=ord)
            {
                TB_inner(ord,trefftzbasis,coeffnum,basis, dim-1, tracker);
                coeffnum(dim-1)++;
            }
        }
        else
        {
            int sum=0;
            for(int i=0;i<D;i++)
                sum += coeffnum(i);
            if(sum<=ord)
            {
                if(tracker >= 0) tracker++;
                int indexmap = IndexMap2(coeffnum, ord);
                if((coeffnum(D-1)==0 || coeffnum(D-1)==1) && tracker>basis)
                {
                    trefftzbasis(basis,indexmap) = 1;
                    tracker = -1;
                }
                else if(coeffnum(D-1)>1)
                {
                    int k = coeffnum(D-1);
                    for(int m=0;m<D-1;m++) //rekursive sum
                    {
                        Vec<D, int> get_coeff = coeffnum;
                        get_coeff[D-1] = get_coeff[D-1] - 2;
                        get_coeff[m] = get_coeff[m] + 2;
                        trefftzbasis( basis, indexmap) += (coeffnum(m)+1) * (coeffnum(m)+2) * trefftzbasis(basis, IndexMap2(get_coeff, ord));
                    }
                    trefftzbasis(basis, indexmap) *= 1.0/(k * (k-1));
                }
            }
        }
    }

    template<int D>
    int TrefftzWaveBasis<D> :: IndexMap2(Vec<D, int> index, int ord)
    {
        int sum=0;
        int temp_size = 0;
        for(int d=0;d<D;d++){
            for(int p=0;p<index(d);p++){
                sum+=BinCoeff(D-1 - d + ord - p - temp_size, ord - p - temp_size);
            }
            temp_size+=index(d);
        }
        return sum;
    }

    template class TrefftzWaveBasis<1>;
    template class TrefftzWaveBasis<2>;
    template class TrefftzWaveBasis<3>;
    template class TrefftzWaveBasis<4>;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef NGS_PYTHON
void ExportTrefftzElement(py::module m)
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
