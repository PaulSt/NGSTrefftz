#ifndef FILE_MONOMIALFESPACE_HPP
#define FILE_MONOMIALFESPACE_HPP

#include "scalarmappedfe.hpp"

namespace ngcomp
{

    class MonomialFESpace : public FESpace
    {
        int D;
        int order;
        size_t ndof;
        int nel;
        int local_ndof;
        int useshift=1;
        shared_ptr<CoefficientFunction> wavespeedcf;
        CSR basismat;

        public:

        MonomialFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

        void SetWavespeed(shared_ptr<CoefficientFunction> awavespeedcf) {
            wavespeedcf=awavespeedcf;
        }

        string GetClassName () const override { return "monomialfespace"; }

        void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

        FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

        size_t GetNDof () const override { return this->ndof; }

        static DocInfo GetDocu ();

        protected:

        template<int D>
        CSR MonomialBasis(int ord) const
        {
            CSR tb;
            const int npoly = BinCoeff(D+1 + ord, ord);
            Matrix<> basis(npoly,npoly);
            basis = 0.0;
            for(int i = 0;i<npoly;i++)
                basis(i,i)=1.0;

            MatToCSR(basis,tb);
            return tb;
        }



        template<int D>
        double Adiam(ElementId ei) const
        {
            LocalHeap lh(1000 * 1000);
            double anisotropicdiam = 0.0;
            auto vertices_index = ma->GetElVertices(ei);

            for(auto vertex1 : vertices_index)
            {
                for(auto vertex2 : vertices_index)
                {
                    Vec<D+1> v1 = ma->GetPoint<D+1>(vertex1);
                    Vec<D+1> v2 = ma->GetPoint<D+1>(vertex2);
                    IntegrationRule ir (ma->GetElType(ei), 0);
                    ElementTransformation & trafo = ma->GetTrafo (ei, lh);
                    MappedIntegrationPoint<D+1,D+1> mip(ir[0], trafo);
                    mip.Point() = v1;
                    double c1 = wavespeedcf->Evaluate(mip);
                    mip.Point() = v2;
                    double c2 = wavespeedcf->Evaluate(mip);

                    anisotropicdiam = max( anisotropicdiam, sqrt( L2Norm2(v1.Range(0,D) - v2.Range(0,D)) + pow(c1*v1(D)-c2*v2(D),2) ) );
                }
            }
            return anisotropicdiam * useshift + (useshift==0);
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
void ExportMonomialFESpace(py::module m);
#endif // NGS_PYTHON

#endif
