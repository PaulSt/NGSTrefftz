#ifndef FILE_QTREFFTZWAVEELEMENT_HPP
#define FILE_QTREFFTZWAVEELEMENT_HPP

#include <fem.hpp>
#include "helpers.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{


    template<int D>
    class QTrefftzWaveBasis{
        public:
            static QTrefftzWaveBasis& getInstance(){
                static QTrefftzWaveBasis ginstance;
                volatile int dummy{};
                return ginstance;
            }

            CSR TB(int ord, FlatMatrix<double> gamma, int basistype = 0);

        private:
            QTrefftzWaveBasis()= default;
            ~QTrefftzWaveBasis()= default;
            QTrefftzWaveBasis(const QTrefftzWaveBasis&)= delete;
            QTrefftzWaveBasis& operator=(const QTrefftzWaveBasis&)= delete;

            mutex gentrefftzbasis;
            std::map<std::string,CSR> gtbstore;
    };


    template <int D>
    class QTrefftzWaveFE : public ScalarMappedElement<D+1>
    {
        private:
            const int ord;
            const int npoly;
            ELEMENT_TYPE eltype;
            int basistype;
            Matrix<double> gamma;

        public:
            QTrefftzWaveFE(Matrix<double> agamma, int aord = 1, Vec<D+1> aelcenter = 0, double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG, int abasistype = 0)
                : ScalarMappedElement<D+1>(BinCoeff(D + aord, aord) + BinCoeff(D + aord-1, aord-1), aord),
                ord(aord),
                npoly(BinCoeff(D+1 + ord, ord)),
                eltype(aeltype),
                basistype(abasistype),
                gamma(agamma)
            {

                this->c = 1.0;
                this->elsize = aelsize;
                this->elcenter = aelcenter;

                static Timer timerbasis("basis",2);
                timerbasis.Start();
                for(int i=0;i<aord-1;i++)
                    for(int j=0;j<aord-1;j++)
                        gamma(i,j) *= pow(aelsize/2.0,i+j);
                this->localmat = QTrefftzWaveBasis<D>::getInstance().TB(ord,gamma);
                timerbasis.Stop();
            }

            double GetWavespeed() const { return gamma(0); }
            //void SetWavespeed(double wavespeed) {gamma(0) = wavespeed;}

            virtual ELEMENT_TYPE ElementType() const { return eltype; }

            using ScalarMappedElement<D+1>::CalcShape;
            using ScalarMappedElement<D+1>::CalcDShape;

            //using ScalarMappedElement<D+1>::CalcMappedDDShape;
            void CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> hddshape) const;
            void CalcDDSpecialShape (const SIMD_BaseMappedIntegrationRule & smir,
                    BareSliceMatrix<SIMD<double>> dshape,
                    BareSliceMatrix<SIMD<double>> wavespeed,
                    BareSliceMatrix<SIMD<double>> mu) const;

            void CalcDDSpecialShape (const SIMD_BaseMappedIntegrationRule & smir,
                    BareSliceMatrix<SIMD<double>> dshape,
                    BareSliceMatrix<SIMD<double>> wavespeed) const
            {
                Matrix<SIMD<double>> mu(1,wavespeed.Dist());
                SIMD<double> a = 1.0;
                mu = a;
                CalcDDSpecialShape (smir,dshape,wavespeed,mu);

            }

            using ScalarMappedElement<D+1>::CalcMappedDShape;
    };



}

#endif // FILE_QTREFFTZWAVEELEMENT_HPP
