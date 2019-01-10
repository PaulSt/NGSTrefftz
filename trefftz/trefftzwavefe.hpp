#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include "helpers.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{
    template <int D>
    class TrefftzWaveFE : public ScalarMappedElement<D>
    {
        private:
            const int ord;
            const int nbasis;
            const int npoly;
            Vec<D> elcenter;
            double elsize;
            float c;
            ELEMENT_TYPE eltype;

        public:
            // TrefftzWaveFE();
            TrefftzWaveFE(int aord = 1, float ac = 1.0, Vec<D> aelcenter = 0, double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG);

            virtual ELEMENT_TYPE ElementType() const { return eltype; }

            using ScalarMappedElement<D>::CalcShape;
            virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;
            virtual void CalcShape (const SIMD_MappedIntegrationRule<D-1,D> & smir, BareSliceMatrix<SIMD<double>> shape) const;

            using ScalarMappedElement<D>::CalcDShape;
            virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;
            virtual void CalcDShape (const SIMD_MappedIntegrationRule<D-1,D> & smir, SliceMatrix<SIMD<double>> dshape) const;

            int GetNBasis() const { return nbasis; }

            //TrefftzWaveFE<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
            //TrefftzWaveFE<D> * SetElSize(double aelsize) {elsize = aelsize; return this;}
            // TrefftzWaveFE<D> * SetWavespeed(float ac) {c = ac; return this;}

        protected:
            void MakeIndices_inner(Matrix<int> &indice, Vec<D, int> numbers, int &count, int ordr, int dim) const;
            Matrix<int> MakeIndices() const;

            constexpr int IndexMap(Vec<D, int> index) const;
            Matrix<double> TrefftzBasis() const;
            Matrix<double> GetDerTrefftzBasis(int der) const;
            Matrix<int> pascal_sym() const;

    };



    class Monomial : public RecursivePolynomial<Monomial>
    {
        public:
            Monomial () { ; }

            template <class S, class T>
                inline Monomial (int n, S x, T && values)
                {
                    Eval (n, x, values);
                }

            template <class S>
                static INLINE double P0(S x)  { return 1.0; }
            template <class S>
                static INLINE S P1(S x)  { return x; }
            template <class S, class Sy>
                static INLINE S P1(S x, Sy y)  { return P1(x); }

            static INLINE double A (int i) { return 1.0; }
            static INLINE double B (int i) { return 0; }
            static INLINE double C (int i) { return 0; }

            static INLINE double CalcA (int i) { return 1.0; }
            static INLINE double CalcB (int i) { return 0; }
            static INLINE double CalcC (int i) { return 0; }

            enum { ZERO_B = 1 };
    };

    static mutex gentrefftzbasis;
    template<int D>
    class TrefftzWaveBasis{
        public:
            static TrefftzWaveBasis& getInstance(){
                static TrefftzWaveBasis instance;
                // volatile int dummy{};
                return instance;
            }

            const Matrix<>* TB(int ord)
            {

                {
                    lock_guard<mutex> lock(gentrefftzbasis);
                    if (tbstore.Size() <= ord)
                    {
                        int oldsize = tbstore.Size();
                        tbstore.SetSize (ord+1);
                        for (int i = oldsize; i <= ord; i++)
                            tbstore[i] = Matrix<>();
                    }

                    if ( tbstore[ord].Height() == 0)
                    {
                        cout << "basis for ord: " << ord << endl;
                        const int nbasis = (BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1));
                        const int npoly = (BinCoeff(D + ord, ord));
                        Matrix<> trefftzbasis(npoly,nbasis);
                        tbstore[ord].SetSize(nbasis,npoly);
                        tbstore[ord] = 0;
                        Vec<D, int>  coeff = 0;
                        int count = 0;
                        for(int b=0;b<nbasis;b++)
                        {
                            int tracker = 0;
                            TB_inner(ord, tbstore[ord], coeff, b, D, tracker);
                        }
                    }

                    if ( tbstore[ord].Height() == 0)
                    {
                        stringstream str;
                        str << "failed to generate trefftz basis of order " << ord << endl;
                        throw Exception (str.str());
                    }
                }
                const Matrix<>* tb =& tbstore[ord];
                return tb;
            }

        private:
            TrefftzWaveBasis()= default;
            ~TrefftzWaveBasis()= default;
            TrefftzWaveBasis(const TrefftzWaveBasis&)= delete;
            TrefftzWaveBasis& operator=(const TrefftzWaveBasis&)= delete;

            Array<Matrix<>> tbstore;

            //once_flag tbonceflag;

            void TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D, int> coeffnum, int basis, int dim, int &tracker)
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
    };

    template class TrefftzWaveBasis<1>;
    template class TrefftzWaveBasis<2>;
    template class TrefftzWaveBasis<3>;
    template class TrefftzWaveBasis<4>;

}


#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement(py::module m);
#endif // NGS_PYTHON


#endif // FILE_TrefftzElement_HPP
