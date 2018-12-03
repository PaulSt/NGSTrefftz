#include "trefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    TrefftzWaveFE<D> ::TrefftzWaveFE(int aord, float ac, ELEMENT_TYPE aeltype, int abasistype)
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

    template<int D>
    void TrefftzWaveFE<D> :: CalcShape (const SIMD_MappedIntegrationRule<D-1,D> & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D,D+1>&> (mir);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<D,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint(D-1) *= c;
            Matrix<SIMD<double>> coeff(TrefftzBasis());

            for(int j=ord;j>0;j--)
                for(int i=0;i<D;i++)
                    for(int k=pascal(i+1,j)-1; k>=0; k--)
                        coeff.Row( pascal(D+1,j-1)+k ) += cpoint[i] * coeff.Row( pascal(D+1,j)+pascal(i,j+1)+k );

            for(int b = 0; b < nbasis; b++) shape.Col(imip)(b) = coeff.Row(0)(b);
        }
    }


    template<int D>
    void TrefftzWaveFE<D> :: CalcDShape (const SIMD_MappedIntegrationRule<D-1,D> & smir,
                                         SliceMatrix<SIMD<double>> dshape) const
    {
        //auto & smir = static_cast<const SIMD_MappedIntegrationRule<D,D+1>&> (mir);
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<D,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[D-1] *= c;

            for(int d=0;d<D;d++){//loop over derivatives/dimensions
                Matrix<SIMD<double>> coeff(GetDerTrefftzBasis(d));
                for(int j=ord-1;j>0;j--)
                    for(int i=0;i<D;i++)
                        for(int k=pascal(i+1,j)-1; k>=0; k--)
                            coeff.Row( pascal(D+1,j-1)+k ) += cpoint[i] * coeff.Row( pascal(D+1,j)+pascal(i,j+1)+k );
                for(int n=0;n<nbasis;n++)
                    dshape(n*D+d,imip) = coeff(0,n) * (d==D-1 ? c : 1);
            }
        }
        dshape *= (2.0/elsize); //inner derivative
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

    template<int D>
    void TrefftzWaveFE<D> :: TB_inner(Matrix<> &trefftzbasis, Vec<D, int> coeffnum, int basis, int ordr, int dim, int &tracker) const
    {
        if (dim>0)
        {
            while(coeffnum(dim-1)<=ordr)
            {
                TB_inner(trefftzbasis,coeffnum,basis, ordr, dim-1, tracker);
                coeffnum(dim-1)++;
            }
        }
        else
        {
            int sum=0;
            for(int i=0;i<D;i++)
                sum += coeffnum(i);
            if(sum<=ordr)
            {
                if(tracker >= 0) tracker++;
                int indexmap = IndexMap(coeffnum);
                if((coeffnum(D-1)==0 || coeffnum(D-1)==1) && tracker>basis)
                {
                    trefftzbasis(indexmap,basis) = 1;
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
                        trefftzbasis( indexmap, basis ) += (coeffnum(m)+1) * (coeffnum(m)+2) * trefftzbasis( IndexMap(get_coeff), basis );
                    }
                    trefftzbasis( indexmap, basis) *= 1.0/(k * (k-1));
                }
            }
        }
    }
    template<int D>
    Matrix<> TrefftzWaveFE<D> :: TB() const
    {
        Matrix<> trefftzbasis(npoly,nbasis);
        trefftzbasis = 0;
        Vec<D, int>  coeff = 0;
        int count = 0;
        for(int b=0;b<nbasis;b++)
        {
            int tracker = 0;
            TB_inner(trefftzbasis, coeff, b, ord, D, tracker);
        }
        return trefftzbasis;
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
