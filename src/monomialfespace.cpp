#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <python_comp.hpp>
#include <fem.hpp>
#include <multigrid.hpp>

#include "monomialfespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{
    MonomialFESpace :: MonomialFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
        : FESpace (ama, flags)
    {
        type="monomialfespace";

        D = ma->GetDimension() - 1;

        order = int(flags.GetNumFlag ("order", 3));
        useshift = flags.GetNumFlag("useshift",1);

        this->local_ndof = BinCoeff(D+1 + order, order);
        this->nel = ma->GetNE();
        this->ndof = local_ndof * nel;

        switch (D)
        {
            case 1:
            {
                evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<2>>>();
                flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<2>>>();
                additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
                basismat = MonomialBasis<1>(order);
                break;
            }
            case 2:
            {
                evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<3>>>();
                flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<3>>>();
                additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
                basismat = MonomialBasis<2>(order);
                break;
            }
        }
    }


    void MonomialFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
        dnums.SetSize(0);
        if (!DefinedOn (ei) || ei.VB() != VOL) return;
        for(size_t j = ei.Nr()*local_ndof; j<local_ndof*(ei.Nr()+1); j++)
            dnums.Append (j);
    }

    template<int D>
    CSR MonomialFESpace :: MonomialBasis(int ord) const
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


    FiniteElement & MonomialFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
        Ngs_Element ngel = ma->GetElement(ei);
        ELEMENT_TYPE eltype = ngel.GetType();

        if (ei.IsVolume())
        {
            switch (ma->GetElType(ei)) {
                case ET_POINT:
                case ET_SEGM:
                {
                    throw Exception("illegal dim for space-time element");
                    break;
                }
                case ET_QUAD:
                case ET_TRIG:
                {
                    return *(new (alloc) ScalarMappedElement<2>(local_ndof,order,basismat,eltype,ElCenter<1>(ei),Adiam<1>(ei)));
                    break;
                }
                case ET_HEX:
                case ET_PRISM:
                case ET_PYRAMID:
                case ET_TET:
                {
                    return *(new (alloc) ScalarMappedElement<3>(local_ndof,order,basismat,eltype,ElCenter<2>(ei),Adiam<2>(ei)));
                    break;
                }
            }
        }
        try
        {
            return SwitchET<ET_POINT,ET_SEGM,ET_TRIG,ET_QUAD>
                (eltype,
                 [&alloc] (auto et) -> FiniteElement&
                 { return * new (alloc) DummyFE<et.ElementType()>; });
        }
        catch (Exception &e)
        {
            throw Exception("illegal element type in Trefftz::GetSurfaceFE");
        }
    }


    template<int D>
    double MonomialFESpace :: Adiam(ElementId ei) const
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
    Vec<D+1> MonomialFESpace :: ElCenter(ElementId ei) const
    {
        Vec<D+1> center = 0;
        auto vertices_index = ma->GetElVertices(ei);
        for(auto vertex : vertices_index) center += ma->GetPoint<D+1>(vertex);
        center *= (1.0/vertices_index.Size()) * useshift;
        return center;
    }

    DocInfo MonomialFESpace :: GetDocu ()
    {
        auto docu = FESpace::GetDocu();
        docu.Arg("useshift") = "bool = True\n"
            "  use shift of basis functins to element center and scale them";
        docu.Arg("gamma") = "bool = True\n"
            "  use shift of basis functins to element center and scale them";
        docu.Arg("basistype") = "bool = True\n"
            "  use shift of basis functins to element center and scale them";
        docu.Arg("wavespeed") = "bool = True\n"
            "  use shift of basis functins to element center and scale them";
        return docu;
    }

    /*
       register fe-spaces
       Object of type MonomialFESpace can be defined in the pde-file via
       "define fespace v -type=Monomialfespace"
       */
    static RegisterFESpace<MonomialFESpace> initi_monomial ("monomialfespace");

}


#ifdef NGS_PYTHON
void ExportMonomialFESpace(py::module m)
{
    using namespace ngcomp;

    ExportFESpace<MonomialFESpace>(m, "monomialfespace")
        .def("GetDocu", &MonomialFESpace::GetDocu)
        .def("GetNDof", &MonomialFESpace::GetNDof)
        .def("SetWavespeed", &MonomialFESpace::SetWavespeed)
        ;
}
#endif // NGS_PYTHON
