#include <comp.hpp> // provides FESpace, ...
#include <python_comp.hpp>

#include "monomialfespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{
  MonomialFESpace ::MonomialFESpace (shared_ptr<MeshAccess> ama,
                                     const Flags &flags, bool)
      : FESpace (ama, flags)
  {
    type = "monomialfespace";

    D = ma->GetDimension ();

    order = int (flags.GetNumFlag ("order", 3));
    useshift = flags.GetNumFlag ("useshift", 1);
    usescale = flags.GetNumFlag ("usescale", 1);

    this->local_ndof = BinCoeff (D + order, order);
    this->nel = ma->GetNE ();
    this->ndof = local_ndof * nel;

    switch (D)
      {
      case 2:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<2>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<2>>> ();
          additional_evaluators.Set (
              "hesse",
              make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
          basismat = MonomialBasis<2> (order);
          break;
        }
      case 3:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<3>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<3>>> ();
          additional_evaluators.Set (
              "hesse",
              make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
          basismat = MonomialBasis<3> (order);
          break;
        }
      }
  }

  void MonomialFESpace ::Update ()
  {
    FESpace::Update ();
    nel = ma->GetNE ();
    ndof = local_ndof * nel;
    SetNDof (ndof);
    UpdateCouplingDofArray ();
  }

  void MonomialFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize (0);
    if (!DefinedOn (ei) || ei.VB () != VOL)
      return;
    for (size_t j = ei.Nr () * local_ndof; j < local_ndof * (ei.Nr () + 1);
         j++)
      dnums.Append (j);
  }

  void MonomialFESpace ::UpdateCouplingDofArray ()
  {
    ctofdof.SetSize (ndof);
    for (auto i : Range (ma->GetNE ()))
      {
        bool definedon = DefinedOn (ElementId (VOL, i));
        Array<DofId> dofs;
        GetDofNrs (i, dofs);
        for (auto r : dofs)
          ctofdof[r] = definedon ? LOCAL_DOF : UNUSED_DOF;
      }
  }

  FiniteElement &MonomialFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement (ei);
    ELEMENT_TYPE eltype = ngel.GetType ();

    if (ei.IsVolume ())
      {
        switch (D)
          {
          case 0:
          case 1:
            {
              throw Exception ("dim not supported in MonomialFESpace");
              break;
            }
          case 2:
            {
              Vec<2> scale = 1.0;
              if (usescale == 2)
                scale = { 1.0 / ElSize<2> (ei, { 1.0, 0 }),
                          1.0 / ElSize<2> (ei, { 0, 1.0 }) };
              else if (usescale != 0)
                scale = 1.0 / ElSize<2> (ei);
              return *(new (alloc) ScalarMappedElement<2> (
                  local_ndof, order, basismat, eltype, ElCenter<2> (ei),
                  scale));
              break;
            }
          case 3:
            {
              Vec<3> scale = 1.0;
              if (usescale == 2)
                {
                  double hx = ElSize<3> (ei, { 1.0, 1.0, 0 });
                  double ht = ElSize<3> (ei, { 0, 0, 1.0 });
                  scale = { 1.0 / hx, 1.0 / hx, 1.0 / ht };
                }
              else if (usescale != 0)
                scale = 1.0 / ElSize<3> (ei);
              return *(new (alloc) ScalarMappedElement<3> (
                  local_ndof, order, basismat, eltype, ElCenter<3> (ei),
                  scale));
              break;
            }
          }
      }
    try
      {
        return SwitchET<ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD> (
            eltype, [&alloc] (auto et) -> FiniteElement & {
              return *new (alloc) DummyFE<et.ElementType ()>;
            });
      }
    catch (Exception &e)
      {
        throw Exception ("illegal element type in Trefftz::GetSurfaceFE");
      }
  }

  DocInfo MonomialFESpace ::GetDocu ()
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
     Object of type MonomialFESpace can be defined in the pde-file via
     "define fespace v -type=Monomialfespace"
     */
  static RegisterFESpace<MonomialFESpace> initi_monomial ("monomialfespace");
}

#ifdef NGS_PYTHON
void ExportMonomialFESpace (py::module m)
{
  using namespace ngcomp;

  ExportFESpace<MonomialFESpace> (m, "monomialfespace")
      .def ("GetDocu", &MonomialFESpace::GetDocu)
      .def ("SetCoeff", &MonomialFESpace::SetCoeff);
}
#endif // NGS_PYTHON
