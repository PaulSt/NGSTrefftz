#include "trefftzwavefe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    TrefftzWaveFE<D> :: TrefftzWaveFE(int aord, float ac, Vec<(D+1)> aelcenter, double aelsize, ELEMENT_TYPE aeltype, int abasistype)
    : ScalarMappedElement<D+1>(BinCoeff(D + aord, aord) + BinCoeff(D + aord-1, aord-1), aord),
    ord(aord),
    npoly(BinCoeff(D+1 + ord, ord)),
    eltype(aeltype),
    basistype(abasistype)
    { 
        this->c = ac;
        this->localmat = TrefftzWaveBasis<D>::getInstance().TB(ord);
        this->elsize = aelsize;
        this->elcenter = aelcenter;
    }


    template class TrefftzWaveFE<1>;
    template class TrefftzWaveFE<2>;
    template class TrefftzWaveFE<3>;
    //template class TrefftzWaveFE<4>;


    template<int D>
    const CSR TrefftzWaveBasis<D> :: TB(int ord)
    {
        return tbstore[ord];
    }

    template<int D>
    void TrefftzWaveBasis<D> :: CreateTB(int ord, int basistype)
    {
        //cout << "creating tb store for order " << ord << endl;
        //if (tbstore.Size() < ord)
        //{
        int oldsize = tbstore.Size();
        tbstore.SetSize (ord+1);
        for (int i = oldsize; i <= ord; i++)
            tbstore[i] = CSR();

        //if ( tbstore[ord][0].Size() == 0)
        //{
        const int ndof = (BinCoeff(D + ord, ord) + BinCoeff(D + ord-1, ord-1));
        const int npoly = (BinCoeff(D+1 + ord, ord));
        Matrix<> trefftzbasis(ndof,npoly);
        trefftzbasis = 0;
        Vec<D+1, int>  coeff = 0;
        int count = 0;
        for(int b=0;b<ndof;b++)
        {
            int tracker = 0;
            TB_inner(ord, trefftzbasis, coeff, b, D+1, tracker, basistype);
        }
        MatToCSR(trefftzbasis,tbstore[ord]);
        //}
        //}

        //if ( tbstore[ord][0].Size() == 0)
        if ( tbstore.Size() < ord)
        {
            stringstream str;
            str << "failed to generate trefftz basis of order " << ord << endl;
            throw Exception (str.str());
        }
    }


    template<int D>
    void TrefftzWaveBasis<D> :: TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D+1, int> coeffnum, int basis, int dim, int &tracker, int basistype, double wavespeed)
    {
        if (dim>0)
        {
            while(coeffnum(dim-1)<=ord)
            {
                TB_inner(ord,trefftzbasis,coeffnum,basis, dim-1, tracker, basistype, wavespeed);
                coeffnum(dim-1)++;
            }
        }
        else
        {
            int sum=0;
            for(int i=0;i<D+1;i++)
                sum += coeffnum(i);
            if(sum<=ord)
            {
                if(tracker >= 0) tracker++;
                int indexmap = IndexMap2(coeffnum, ord);
                int k = coeffnum(D);
                if(k==0 || k==1)
                {
                    switch (basistype) {
                        case 0:
                            if(tracker>basis)
                            {
                                //trefftzbasis( i, setbasis++ ) = 1.0; //set the l-th coeff to 1
                                trefftzbasis(basis,indexmap) = 1;
                                tracker = -1;
                            }
                            //i += ndof-1;	//jump to time = 2 if i=0
                            break;
                        case 1:
                            if((k == 0 && basis < BinCoeff(D + ord, ord)) || (k == 1 && basis >= BinCoeff(D + ord, ord))){
                                trefftzbasis( basis,indexmap ) = 1;
                                for(int exponent :  coeffnum.Range(0,D)) trefftzbasis( basis,indexmap ) *= LegCoeffMonBasis(basis,exponent);}
                            break;
                        case 2:
                            if((k == 0 && basis < BinCoeff(D + ord, ord)) || (k == 1 && basis >= BinCoeff(D + ord, ord))){
                                trefftzbasis( basis,indexmap ) = 1;
                                for(int exponent :  coeffnum.Range(0,D)) trefftzbasis( basis,indexmap ) *= ChebCoeffMonBasis(basis,exponent);}
                            break;
                    }
                }
                else if(coeffnum(D)>1)
                {
                    for(int m=0;m<D;m++) //rekursive sum
                    {
                        Vec<D+1, int> get_coeff = coeffnum;
                        get_coeff[D] = get_coeff[D] - 2;
                        get_coeff[m] = get_coeff[m] + 2;
                        trefftzbasis( basis, indexmap) += (coeffnum(m)+1) * (coeffnum(m)+2) * trefftzbasis(basis, IndexMap2(get_coeff, ord));
                    }
                    trefftzbasis(basis, indexmap) *= wavespeed*wavespeed/(k * (k-1));
                }
            }
        }
    }

    template<int D>
    int TrefftzWaveBasis<D> :: IndexMap2(Vec<D+1, int> index, int ord)
    {
        int sum=0;
        int temp_size = 0;
        for(int d=0;d<D+1;d++){
            for(int p=0;p<index(d);p++){
                sum+=BinCoeff(D - d + ord - p - temp_size, ord - p - temp_size);
            }
            temp_size+=index(d);
        }
        return sum;
    }

    template class TrefftzWaveBasis<1>;
    template class TrefftzWaveBasis<2>;
    template class TrefftzWaveBasis<3>;
    //template class TrefftzWaveBasis<4>;
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
