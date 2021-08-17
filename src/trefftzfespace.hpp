#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP

#include "scalarmappedfe.hpp"
#include <boost/math/special_functions/airy.hpp>

namespace ngcomp
{

    struct GenericSqrt {
        template <typename T> T operator() (T x) const { return sqrt(x); }
        static string Name() { return "sqrt"; }
        void DoArchive(Archive& ar) {}
    };

    class PolBasis {};

    template<int D>
    class TWaveBasis : public PolBasis
    {
        static void TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D+1, int> coeffnum, int basis, int dim, int &tracker, int basistype, double wavespeed = 1.0);
        public:
        TWaveBasis() {;}
        static CSR Basis(int ord, int basistype = 0);
        static int IndexMap2(Vec<D+1, int> index, int ord);
    };

    template<int D>
    class QTWaveBasis : public PolBasis
    {
        mutex gentrefftzbasis;
        std::map<string,CSR> gtbstore;
        public:
        QTWaveBasis() {;}
        CSR Basis(int ord, Vec<D+1> ElCenter, Matrix<shared_ptr<CoefficientFunction>> GGder, Matrix<shared_ptr<CoefficientFunction>> BBder, double elsize = 1.0, int basistype=0);
    };

    class TrefftzFESpace : public FESpace
    {
        int D;
        int order;
        size_t ndof;
        int nel;
        int nvert;
        int local_ndof;
        float c=1;
        int useshift=1;
        int usescale=1;
        int useqt = 0;
        int basistype;
        shared_ptr<CoefficientFunction> wavespeedcf=nullptr;
        Matrix<shared_ptr<CoefficientFunction>> GGder;
        Matrix<shared_ptr<CoefficientFunction>> BBder;
        CSR basismat;
        PolBasis* basis;

        public:
            TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);
            void SetWavespeed(shared_ptr<CoefficientFunction> awavespeedcf, shared_ptr<CoefficientFunction> aBBcf = nullptr);
            string GetClassName () const override { return "trefftz"; }
            void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
            FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
            size_t GetNDof () const override { return ndof; }
            static DocInfo GetDocu ();

        protected:
            template<int D>
            double Adiam(ElementId ei, double c) const;
            template<int D>
            double Adiam(ElementId ei, shared_ptr<CoefficientFunction> c) const;
            template<int D>
            Vec<D+1> ElCenter(ElementId ei) const;
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

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace(py::module m);
#endif // NGS_PYTHON

#endif
