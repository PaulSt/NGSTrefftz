#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP

#include "scalarmappedfe.hpp"

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
        static CSR Basis(int ord, int basistype = 0, int fosystem = 0);
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
            void SetWavespeed(shared_ptr<CoefficientFunction> awavespeedcf, shared_ptr<CoefficientFunction> aBBcf = nullptr, shared_ptr<CoefficientFunction> aGGcf = nullptr);
            string GetClassName () const override { return "trefftz"; }
            void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
            FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
            size_t GetNDof () const override { return ndof; }
            static DocInfo GetDocu ();

        protected:
            template<int D>
            double Adiam(ElementId ei, double c) const
        {
            double anisotropicdiam = 0.0;
            auto vertices_index = ma->GetElVertices(ei);
            for(auto vertex1 : vertices_index)
            {
                for(auto vertex2 : vertices_index)
                {
                    Vec<D+1> v1 = ma->GetPoint<D+1>(vertex1);
                    Vec<D+1> v2 = ma->GetPoint<D+1>(vertex2);
                    anisotropicdiam = max( anisotropicdiam, sqrt( L2Norm2(v1.Range(0,D) - v2.Range(0,D)) + pow(c*(v1(D)-v2(D)),2) ) );
                }
            }
            return anisotropicdiam * usescale + (usescale==0);
        }
            template<int D>
            Vec<D+1> ElCenter(ElementId ei) const
        {
            Vec<D+1> center = 0;
            auto vertices_index = ma->GetElVertices(ei);
            for(auto vertex : vertices_index) center += ma->GetPoint<D+1>(vertex);
            center *= (1.0/vertices_index.Size()) * useshift;
            return center;
        }

    };



}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace(py::module m);
#endif // NGS_PYTHON

#endif
