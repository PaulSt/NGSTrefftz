#include <comp.hpp> // provides FESpace, ...
#include <python_comp.hpp>

#include "pufe.hpp"
#include "pufespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{
  PUFESpace ::PUFESpace (shared_ptr<MeshAccess> ama, const Flags &flags, bool)
      : FESpace (ama, flags)
  {
    type = "pufespace";

    D = ma->GetDimension ();

    order = int (flags.GetNumFlag ("order", 3));
    useshift = flags.GetNumFlag ("useshift", 1);
    usescale = flags.GetNumFlag ("usescale", 1);

    switch (D)
      {
      case 1:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<1>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<1>>> ();
          basismat = MonomialBasis<1> (order);
          break;
        }
      case 2:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<2>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<2>>> ();
          basismat = MonomialBasis<2> (order);
          break;
        }
      case 3:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<3>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<3>>> ();
          basismat = MonomialBasis<3> (order);
          break;
        }
      }
  }

  void PUFESpace ::Update ()
  {
    this->local_ndof = BinCoeff (D + order, order);
    int nvert = ma->GetNV ();

    size_t ndof = local_ndof * nvert;
    SetNDof (ndof);

    // FESpace::Update ();
    UpdateCouplingDofArray ();
  }

  void PUFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize0 ();
    if (!DefinedOn (ei)) // || ei.VB () != VOL)
      return;
    // Ngs_Element ngel = ma->GetElement(ei);
    // dnums = ngel.Vertices();
    for (auto v : ma->GetElVertices (ei))
      for (size_t j = 0; j < local_ndof; j++)
        dnums.Append (v * local_ndof + j);
  }

  void PUFESpace ::UpdateCouplingDofArray ()
  {
    ctofdof.SetSize (ndof);
    for (auto i : Range (ma->GetNE ()))
      {
        bool definedon = DefinedOn (ElementId (VOL, i));
        Array<DofId> dofs;
        GetDofNrs (i, dofs);
        for (auto r : dofs)
          ctofdof[r] = definedon ? WIREBASKET_DOF : UNUSED_DOF;
      }
  }

  FiniteElement &PUFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement (ei);
    ELEMENT_TYPE eltype = ngel.GetType ();

    if (ei.IsVolume ())
      {
        switch (ma->GetElType (ei))
          {
          case ET_SEGM:
            {
              return *(new (alloc) PUFElement<1> (
                  local_ndof * 2, order, basismat, eltype, ElVertices<1> (ei),
                  Adiam<1> (ei)));
              break;
            }
          case ET_TRIG:
            {
              return *(new (alloc) PUFElement<2> (
                  local_ndof * 3, order, basismat, eltype, ElVertices<2> (ei),
                  Adiam<2> (ei)));
              break;
            }
          case ET_TET:
            {
              return *(new (alloc) PUFElement<3> (
                  local_ndof * 4, order, basismat, eltype, ElVertices<3> (ei),
                  Adiam<3> (ei)));
              break;
            }
          default:
            throw Exception ("illegal el for pufespace");
          }
      }
    else
      throw Exception ("Boundary elements not implemented yet!");
    // try
    //{
    // return SwitchET<ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD> (
    // eltype, [&alloc] (auto et) -> FiniteElement & {
    // return *new (alloc) DummyFE<et.ElementType ()>;
    //});
    //}
    // catch (Exception &e)
    //{
    // throw Exception ("illegal element type in Trefftz::GetSurfaceFE");
    //}
  }

  DocInfo PUFESpace ::GetDocu ()
  {
    auto docu = FESpace::GetDocu ();
    docu.Arg ("useshift") = "bool = True\n"
                            "  shift of basis functins to element center";
    docu.Arg ("usescale") = "bool = True\n"
                            "  scale element basis functions with diam";
    return docu;
  }

  /*
     register fe-spaces
     Object of type PUFESpace can be defined in the pde-file via
     "define fespace v -type=PUFESpace"
     */
  static RegisterFESpace<PUFESpace> initi_pufem ("pufespace");
}

#ifdef NGS_PYTHON
void ExportPUFESpace (py::module m)
{
  using namespace ngcomp;

  ExportFESpace<PUFESpace> (m, "PUFESpace")
      .def ("GetDocu", &PUFESpace::GetDocu);
}
#endif // NGS_PYTHON
