#include "trefftzwavebasis.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    void TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D, int> coeffnum, int basis, int dim, int &tracker) 
    {
        if (dim>0)
        {
            while(coeffnum(dim-1)<=ord)
            {
                TB_inner<D>(ord,trefftzbasis,coeffnum,basis, dim-1, tracker);
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
                int indexmap = IndexMap2<D>(coeffnum, ord);
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
                        trefftzbasis( basis, indexmap) += (coeffnum(m)+1) * (coeffnum(m)+2) * trefftzbasis(basis, IndexMap2<D>(get_coeff, ord));
                    }
                    trefftzbasis(basis, indexmap) *= 1.0/(k * (k-1));
                }
            }
        }
    }

    template<int D>
    const Matrix<>*  TB(int ord) 
    {
        //static Matrix<> *tbstore[15] = {NULL};
        static Matrix<> trefftzbasis;
        std::mutex m;
        std::lock_guard<std::mutex> lock{m};

        //if(tbstore[ord]==NULL)
        if (trefftzbasis.Height()==0)
        {
            const int nbasis = (BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1));
            const int npoly = (BinCoeff(D + ord, ord));
            //Matrix<> trefftzbasis(npoly,nbasis);
            trefftzbasis.SetSize(nbasis,npoly);
            trefftzbasis = 0;
            Vec<D, int>  coeff = 0;
            int count = 0;
            for(int b=0;b<nbasis;b++)
            {
                int tracker = 0;
                TB_inner<D>(ord, trefftzbasis, coeff, b, D, tracker);
            }
            //*tbstore[ord] = trefftzbasis;
        }
        //return tbstore[ord];
        const Matrix<>* tb =&trefftzbasis;
        return tb;
    }

    template<int D>
    int IndexMap2(Vec<D, int> index, int ord) 
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

    template const Matrix<>* TB<2>(int ord) ;
}
