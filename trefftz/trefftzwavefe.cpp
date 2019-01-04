#include "trefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"
#include "trefftzwavebasis.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    TrefftzWaveFE<D> :: TrefftzWaveFE(int aord, float ac, ELEMENT_TYPE aeltype, int abasistype)
    : ScalarMappedElement<D>(BinCoeff(D-1 + aord, aord) + BinCoeff(D-1 + aord-1, aord-1), aord),
    ord(aord),
    c(ac),
    nbasis(BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1)),
    npoly(BinCoeff(D + ord, ord)),
    basistype(abasistype),
    eltype(aeltype),
    pascal(pascal_sym())
    {;}

    template<int D>
    void TrefftzWaveFE<D> :: CalcShape (const SIMD_MappedIntegrationRule<D-1,D> & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (mir);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<D,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;
            SIMD<double> x = cpoint[0], t = cpoint[1];

            STACK_ARRAY(SIMD<double>, mem, D*(ord+1));
            Vec<D,SIMD<double>*> polxt;
            for(size_t d=0;d<D;d++)
            {
                polxt[d] = &mem[d*(ord+1)];
                Monomial (ord, cpoint[d], polxt[d]);
            }

            Matrix<double> localmat = *TB<2>(ord);
            Vector<SIMD<double>> tempshape(nbasis);
            Vector<SIMD<double>> pol(npoly);
            pol= 1;
            //tempshape=0;
            //
            Vec<D+1,int> ind = 0; // if "n" is not known before hand, then this array will need to be created dynamically.
            int p = 0; //Used to increment all of the indicies correctly, at the end of each loop.
            int ii = 0; // keep track of runs
            while (ind[D]==0) 
            {
                int sum=0;
                for(int i=0;i<=D;i++)
                    sum += ind(i);
                if(sum<=ord)
                {
                    for(int d=0;d<D;d++)
                        pol[ii] = pol[ii] * polxt[d][ind[d]];
                    ii++;
                    //cout << ind << endl;
                }
                // increment all of the indicies correctly.
                ind[0]++;
                p=0;
                while(ind[p]>ord) {
                    ind[p]=0;
                    ind[++p]++; //increase p by 1, and increase the next (p+1)th index
                    //if(ind[p] != ord)
                    //p=0; // only reset p when it actually needs to be reset!
                }
            }
            //cout << endl << endl;
            //for (size_t i = 0, ii = 0; i <= ord; i++)
            //for (size_t j = 0; j <= ord-i; j++)
            //{
            //pol[ii++] = polx[i] * polt[j];
            ////tempshape += pol * localmat.Col(ii++);
            //}
            tempshape = localmat * pol;
            for(int b=0;b<nbasis;b++) shape.Col(imip)(b) = tempshape(b);
        }
    }

    template<>
    void TrefftzWaveFE<2> :: CalcDShape (const SIMD_MappedIntegrationRule<1,2> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;
            SIMD<double> x = cpoint[0], t = cpoint[1];

            STACK_ARRAY(SIMD<double>, mem, 2*ord+2);
            SIMD<double> * polx = &mem[0];
            SIMD<double> * polt = &mem[ord+1];

            Monomial (ord, x, polx);
            Monomial (ord, t, polt);

            Matrix<> localmat = *TB<2>(ord);
            Vector<SIMD<double>> tempdshape(nbasis);
            for(int d=0;d<2;d++)
            {
                tempdshape=0;
                for (size_t i = 0, ii = 0; i <=ord; i++)
                    for (size_t j = 0; j <= ord-i; j++)
                    {
                        ii++;
                        if((d==0&&j==0)||(d==1&&i==0)) continue;
                        SIMD<double> pol = polx[j-(d==0)] * polt[i-(d==1)] * (d==0?j:i);
                        tempdshape += pol * localmat.Col(ii-1);
                    }
                //dshape.Col(d) = tempdshape;
                for(int n=0;n<nbasis;n++)
                    dshape(n*2+d,imip) = tempdshape(n) * (d==1 ? c : 1);
            }
        }
        //dshape.Col(1) *= c; //inner derivative
        dshape *= (2.0/elsize); //inner derivative
    }

    template<int D>
    void TrefftzWaveFE<D> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
        Vec<D> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[D-1] *= c;
        Matrix<> coeff(TrefftzBasis());

        for(int j=ord;j>0;j--)
            for(int i=0;i<D;i++)
                for(int k=pascal(i+1,j)-1; k>=0; k--)
                    coeff.Row( pascal(D+1,j-1)+k ) += cpoint[i] * coeff.Row( pascal(D+1,j)+pascal(i,j+1)+k );

        for(int b = 0; b < nbasis; b++) shape(b) = coeff.Row(0)(b);
    }


    template<int D>
    void TrefftzWaveFE<D> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {
        Vec<D> cpoint = mip.GetPoint();
        cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[D-1] *= c;

        for(int d=0;d<D;d++){//loop over derivatives/dimensions
            Matrix<double> coeff(GetDerTrefftzBasis(d));
            for(int j=ord-1;j>0;j--)
                for(int i=0;i<D;i++)
                    for(int k=pascal(i+1,j)-1; k>=0; k--)
                        coeff.Row( pascal(D+1,j-1)+k ) += cpoint[i] * coeff.Row( pascal(D+1,j)+pascal(i,j+1)+k );
            dshape.Col(d) = coeff.Row(0);
        }

        dshape.Col(D-1) *= c; //inner derivative
        dshape *= (2.0/elsize); //inner derivative
    }

    //template<int D>
    //void TrefftzWaveFE<D> :: CalcShape (const SIMD_MappedIntegrationRule<D-1,D> & smir,
    //BareSliceMatrix<SIMD<double>> shape) const
    //{
    ////auto & smir = static_cast<const SIMD_MappedIntegrationRule<D,D+1>&> (mir);
    //for (int imip = 0; imip < smir.Size(); imip++)
    //{
    //Vec<D,SIMD<double>> cpoint = smir[imip].GetPoint();
    //cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint(D-1) *= c;
    //Matrix<SIMD<double>> coeff(TrefftzBasis());

    //for(int j=ord;j>0;j--)
    //for(int i=0;i<D;i++)
    //for(int k=pascal(i+1,j)-1; k>=0; k--)
    //coeff.Row( pascal(D+1,j-1)+k ) += cpoint[i] * coeff.Row( pascal(D+1,j)+pascal(i,j+1)+k );

    //for(int b = 0; b < nbasis; b++) shape.Col(imip)(b) = coeff.Row(0)(b);
    //}
    //}


    template<int D>
    void TrefftzWaveFE<D> :: CalcDShape (const SIMD_MappedIntegrationRule<D-1,D> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D,D+1>&> (mir);

        constexpr int sw = SIMD<double>::Size();
        static Timer calcdshape("calcdshape",2);
        RegionTimer T(calcdshape);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<D,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[D-1] *= c;

            for(int d=0;d<D;d++){//loop over derivatives/dimensions
                Matrix<SIMD<double>> coeff(GetDerTrefftzBasis(d));
                for(int j=ord-1;j>0;j--)
                    for(int i=0;i<D;i++)
                        for(int k=pascal(i+1,j)-1; k>=0; k--)
                        {
                            coeff.Row( pascal(D+1,j-1)+k ) += cpoint[i] * coeff.Row( pascal(D+1,j)+pascal(i,j+1)+k );
                            //calcdshape.AddFlops(coeff.Width()*sw);
                        }
                for(int n=0;n<nbasis;n++)
                    dshape(n*D+d,imip) = coeff(0,n) * (d==D-1 ? c : 1);
                //calcdshape.AddFlops(nbasis*sw);
            }
        }
        dshape *= (2.0/elsize); //inner derivative
        //calcdshape.AddFlops(nbasis*(D+1)*sw);
    }

    template<int D>
    Matrix<double> TrefftzWaveFE<D> :: GetDerTrefftzBasis(int der) const
    {
        static int order;
        static int btype;
        static Vec<D,Matrix<double>> basisstorage;
        if(order != ord || btype != basistype){
            Matrix<int> indices = MakeIndices();
            for(int d=0;d<D;d++){
                basisstorage[d].SetSize(pascal(D+1,ord), nbasis);
                for(int i=0, count=0;i<npoly;i++)
                {
                    if(indices(i,d) != 0)
                        basisstorage[d].Row(count++) = indices(i,d)*TrefftzBasis().Row(i);
                }
            }
            order = ord;
            btype = basistype;
        }
        return basisstorage[der];
    }

    template<int D>
    Matrix<double> TrefftzWaveFE<D> :: TrefftzBasis() const
    {
        static int order;
        static int btype;
        static Matrix<double> basisstorage;

        if(order != ord || btype != basistype)
        {
            basisstorage.SetSize(npoly,nbasis);

            basisstorage = 0;
            int setbasis = 0;
            Matrix<int> indices = MakeIndices();
            for(int l=0;l<nbasis;l++) //loop over basis functions
            {
                for(int i=0;i<npoly;i++)//loop over indices BinCoeff(D + ord, ord)
                {
                    int k = indices(i,D-1);
                    if(k > 1)
                    {
                        for(int m=0;m<D-1;m++) //rekursive sum
                        {
                            Vec<D, int> get_coeff = indices.Row(i);
                            get_coeff[D-1] = get_coeff[D-1] - 2;
                            get_coeff[m] = get_coeff[m] + 2;
                            basisstorage( i, l ) += (indices(i,m)+1) * (indices(i,m)+2) * basisstorage( IndexMap(get_coeff), l );
                        }
                        basisstorage( i, l ) *= 1.0/(k * (k-1));
                    }
                    else if(k<=1) //time=0 and =1
                    {
                        switch (basistype) {
                            case 0:
                                if(l==0)
                                    basisstorage( i, setbasis++ ) = 1.0; //set the l-th coeff to 1
                                //i += nbasis-1;	//jump to time = 2 if i=0
                                break;
                            case 1:
                                if((k == 0 && l < BinCoeff(D-1 + ord, ord)) || (k == 1 && l >= BinCoeff(D-1 + ord, ord))){
                                    basisstorage( i, l ) = 1;
                                    for(int exponent :  indices.Row(i).Range(0,D-1)) basisstorage( i, l ) *= LegCoeffMonBasis(l,exponent);}
                                break;
                            case 2:
                                if((k == 0 && l < BinCoeff(D-1 + ord, ord)) || (k == 1 && l >= BinCoeff(D-1 + ord, ord))){
                                    basisstorage( i, l ) = 1;
                                    for(int exponent :  indices.Row(i).Range(0,D-1)) basisstorage( i, l ) *= ChebCoeffMonBasis(l,exponent);}
                                break;
                        }
                    }
                }
            }
            order = ord;
            btype = basistype;
        }
        return basisstorage;
    }


    template<int D>
    void TrefftzWaveFE<D> :: MakeIndices_inner(Matrix<int> &indice, Vec<D, int> numbers, int &count, int ordr, int dim) const
    {
        if (dim>0)
        {
            for(int i=0;i<=ordr;i++)
            {
                numbers(dim-1)=i;
                MakeIndices_inner(indice,numbers,count, ordr, dim-1);
            }
        }
        else
        {
            int sum=0;
            for(int i=0;i<D;i++)
                sum += numbers(i);
            if(sum==ordr)
                indice.Row(count++) = numbers;
        }
    }


    template<int D>
    Matrix<int> TrefftzWaveFE<D> :: MakeIndices() const
    {
        Matrix<int> indice(npoly,D);
        Vec<D, int>  numbers = 0;
        int count = 0;
        for(int o=0;o<=ord;o++)
            MakeIndices_inner(indice, numbers, count, o, D);
        return indice;
    }


    template<int D>
    constexpr int TrefftzWaveFE<D> :: IndexMap(Vec<D, int> index) const
    {
        int sum = 0;
        int indexleng = 0;
        for(int r=0; r<D; r++)
        {
            indexleng += index(r);
            for(int i=0;i<index(r);i++){
                sum+=BinCoeff(indexleng-i+r-1,indexleng-i);}
        }
        sum += BinCoeff(indexleng-1+D,indexleng-1);
        return sum;
    }


    template<int D>
    Matrix<int> TrefftzWaveFE<D> :: pascal_sym() const
    {
        static int order;
        static Matrix<int> pascalstorage;

        if(order != ord)
        {
            pascalstorage.SetSize(D+2,ord+2);
            for (int i = 0; i <= D+1; ++i)
                for (int j = 0; j <= ord+1; ++j)
                    if (i == 0 || j == 0)
                        pascalstorage(i,j) = 0;
                    else if (i == 1 || j == 1)
                        pascalstorage(i,j) = 1;
                    else
                        pascalstorage(i,j) = pascalstorage(i-1,j) + pascalstorage(i,j-1);

            order = ord;
        }
        return pascalstorage;
    }



    template class TrefftzWaveFE<1>;
    template class TrefftzWaveFE<2>;
    template class TrefftzWaveFE<3>;
    template class TrefftzWaveFE<4>;
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
