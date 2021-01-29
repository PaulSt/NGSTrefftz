#ifndef FILE_MONOMIALELEMENT_HPP
#define FILE_MONOMIALELEMENT_HPP

#include <fem.hpp>
#include "helpers.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{


    template<int D>
    class MonomialBasis{
        public:
            MonomialBasis(int ord)
            {
                const int npoly = BinCoeff(D+1 + ord, ord);
                Matrix<> basis(npoly,npoly);
                basis = 0.0;
                for(int i = 0;i<npoly;i++)
                    basis(i,i)=1.0;

                MatToCSR(basis,tb);
            }

            CSR TB() const {return tb;}

        private:
            CSR tb;
    };

    template class MonomialBasis<1>;
    template class MonomialBasis<2>;


    template <int D>
    class MonomialFE : public ScalarMappedElement<D+1>
    {
        private:
            const int ord;
            const int npoly;
            ELEMENT_TYPE eltype;
            Matrix<double> gamma;

        public:
            MonomialFE(int aord = 1, Vec<D+1> aelcenter = 0, double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG)
            : ScalarMappedElement<D+1>( BinCoeff(D+1 + aord, aord), aord),
            ord(aord),
            npoly(BinCoeff(D+1 + aord, aord)),
            eltype(aeltype)
            {
                MonomialBasis<D> Basis(aord);
                this->localmat = Basis.TB();
                this->elsize = aelsize;
                this->elcenter = aelcenter;
                this->c = 1.0;
            }

            virtual ELEMENT_TYPE ElementType() const { return eltype; }

            using ScalarMappedElement<D+1>::CalcShape;
            using ScalarMappedElement<D+1>::CalcDShape;

            using ScalarMappedElement<D+1>::CalcMappedDShape;
            using ScalarMappedElement<D+1>::CalcMappedDDShape;

    };

    template class MonomialFE<1>;
    template class MonomialFE<2>;

}

#endif // FILE_TrefftzGPPWElement_HPP
