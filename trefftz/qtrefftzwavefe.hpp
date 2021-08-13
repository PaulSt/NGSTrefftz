//////////////////////////////////////////
// Quasi-Trefftz space for the equation
//    div(B grad(u)) - G d_tt(u)= 0
// where B=B(x) and G=G(x) are smooth
//////////////////////////////////////////

#ifndef FILE_QTREFFTZWAVEELEMENT_HPP
#define FILE_QTREFFTZWAVEELEMENT_HPP

#include <fem.hpp>
#include "scalarmappedfe.hpp"
#include <cmath>
#include <bla.hpp>
#include <fem.hpp>
#include <boost/math/special_functions/airy.hpp>

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

            CSR TB(int ord, Vec<D+1> ElCenter, Matrix<shared_ptr<CoefficientFunction>> GGder, Matrix<shared_ptr<CoefficientFunction>> BBder, double elsize = 1.0, int basistype=0);

            void Clear(){gtbstore.clear();}

        private:
            QTrefftzWaveBasis()= default;
            ~QTrefftzWaveBasis()= default;
            QTrefftzWaveBasis(const QTrefftzWaveBasis&)= delete;
            QTrefftzWaveBasis& operator=(const QTrefftzWaveBasis&)= delete;

            mutex gentrefftzbasis;
            std::map<string,CSR> gtbstore;
    };


    template <int D>
    class QTrefftzWaveFE : public ScalarMappedElement<D+1>
    {
        private:
            const int ord;
            const int npoly;
            ELEMENT_TYPE eltype;
            int basistype;

        public:
            QTrefftzWaveFE(Matrix<shared_ptr<CoefficientFunction>> aGGder, Matrix<shared_ptr<CoefficientFunction>> aBBder, int aord = 1, Vec<D+1> aelcenter = 0, double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG, int abasistype = 0)
                : ScalarMappedElement<D+1>(BinCoeff(D + aord, aord) + BinCoeff(D + aord-1, aord-1), aord),
                ord(aord),
                npoly(BinCoeff(D+1 + ord, ord)),
                eltype(aeltype),
                basistype(abasistype)
            {

                this->c = 1.0;
                this->elsize = aelsize/2.0;
                this->elcenter = aelcenter;

                static Timer timerbasis("quasiTrefftzbasis");
                timerbasis.Start();
                this->localmat = QTrefftzWaveBasis<D>::getInstance().TB(ord,aelcenter,aGGder,aBBder,aelsize/2.0);
                timerbasis.Stop();
            }

            double GetWavespeed() const { return 0; }

            ELEMENT_TYPE ElementType() const { return eltype; }

            using ScalarMappedElement<D+1>::CalcShape;
            using ScalarMappedElement<D+1>::CalcDShape;


            using ScalarMappedElement<D+1>::CalcMappedDShape;
    };



    struct GenericAiry {
        double operator() (double x) const { return boost::math::airy_ai(x); }
        SIMD<double> operator() (SIMD<double> x) const
        {
            return SIMD<double>([&](int i)->double
                                {
                                    return boost::math::airy_ai(x[i]);
                                });
        }
        template <typename T> T operator() (T x) const
        { throw Exception (string("airy not available for type ")+typeid(T).name());  }
        template<typename T>
        AutoDiff<1,T> operator() (AutoDiff<1,T> x) const { throw Exception ("airy(..) is not complex differentiable"); }
        template<typename T>
        AutoDiffDiff<1,T> operator() (AutoDiffDiff<1,T> x) const { throw Exception ("airy(..) is not complex differentiable"); }
        static string Name() { return "airy"; }
        void DoArchive(Archive& ar) {}
    };

    struct GenericAiryP {
        double operator() (double x) const { return boost::math::airy_ai_prime(x); }
        SIMD<double> operator() (SIMD<double> x) const
        {
            return SIMD<double>([&](int i)->double
                                {
                                    return boost::math::airy_ai_prime(x[i]);
                                });
        }
        template <typename T> T operator() (T x) const
        { throw Exception (string("airy prime not available for type ")+typeid(T).name());  }
        template<typename T>
        AutoDiff<1,T> operator() (AutoDiff<1,T> x) const { throw Exception ("airyp(..) is not complex differentiable"); }
        template<typename T>
        AutoDiffDiff<1,T> operator() (AutoDiffDiff<1,T> x) const { throw Exception ("airyp(..) is not complex differentiable"); }
        static string Name() { return "airyp"; }
        void DoArchive(Archive& ar) {}
    };

}

#endif // FILE_QTREFFTZWAVEELEMENT_HPP
