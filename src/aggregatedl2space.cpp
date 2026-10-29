#include "aggregatedl2space.hpp"

#include <python_comp.hpp>

namespace ngcomp
{
  AggregatedL2Space::AggregatedL2Space (
      shared_ptr<MeshAccess> ama, const Flags &flags,
      shared_ptr<ElementAggregation> aaggregation)
      : MonomialFESpace (ama, Flags (flags).SetFlag (
                                  "order", flags.GetNumFlag ("order", 0))),
        aggregation (aaggregation)
  {
    type = "aggregatedl2";
  }

  void AggregatedL2Space::Update ()
  {
    FESpace::Update ();

    const size_t ne = ma->GetNE (VOL);
    const size_t npatches = aggregation->GetNPatches (true);
    element_to_patch.SetSize (ne);
    element_to_patch = -1;
    centers2.SetSize (npatches);
    centers3.SetSize (npatches);
    scales.SetSize (npatches);

    Array<size_t> els;
    for (size_t patch = 0; patch < npatches; patch++)
      {
        aggregation->GetPatch (patch, els);
        for (size_t el : els)
          element_to_patch[el] = patch;
        if (D == 2)
          UpdateGeometry<2> (patch, els);
        else
          UpdateGeometry<3> (patch, els);
      }

    SetNDof (npatches * local_ndof);
    UpdateCouplingDofArray ();
  }

  void AggregatedL2Space::UpdateCouplingDofArray ()
  {
    ctofdof.SetSize (GetNDof ());
    ctofdof = UNUSED_DOF;
    for (size_t el : Range (ma->GetNE (VOL)))
      if (element_to_patch[el] >= 0 && DefinedOn (ElementId (VOL, el)))
        for (size_t j : Range (local_ndof))
          ctofdof[element_to_patch[el] * local_ndof + j]
              = aggregation->GetElsInNontrivialPatch ()->Test (el)
                    ? WIREBASKET_DOF
                    : LOCAL_DOF;
  }

  void AggregatedL2Space::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize0 ();
    if (ei.VB () != VOL || !DefinedOn (ei) || element_to_patch[ei.Nr ()] < 0)
      return;
    const int first = element_to_patch[ei.Nr ()] * local_ndof;
    for (size_t j : Range (local_ndof))
      dnums.Append (first + j);
  }

  DofPatches AggregatedL2Space::GetDofPatches () const
  {
    size_t npatches = aggregation->GetNPatches (true);
    TableCreator<size_t> elements (npatches);
    Array<size_t> els;
    for (; !elements.Done (); elements++)
      for (size_t p : Range (npatches))
        {
          aggregation->GetPatch (p, els);
          for (size_t el : els)
            elements.Add (p, el);
        }
    DofPatches patches{ elements.MoveTable (),
                        Table<DofId> (npatches, local_ndof) };
    auto dofs = patches.dofs.AsArray ();
    for (size_t d : Range (dofs))
      dofs[d] = d;
    return patches;
  }

  FiniteElement &
  AggregatedL2Space::GetFE (ElementId ei, Allocator &alloc) const
  {
    if (ei.VB () != VOL)
      return MonomialFESpace::GetFE (ei, alloc);
    const int patch = element_to_patch[ei.Nr ()];
    const ELEMENT_TYPE et = ma->GetElType (ei);
    if (patch < 0 || !DefinedOn (ei))
      return SwitchET<ET_TRIG, ET_QUAD, ET_TET, ET_HEX, ET_PRISM, ET_PYRAMID> (
          et, [&alloc] (auto type) -> FiniteElement & {
            return *new (alloc) DummyFE<type.ElementType ()>;
          });
    if (D == 2)
      return *(new (alloc) ScalarMappedElement<2> (local_ndof, order, basismat,
                                                   et, centers2[patch],
                                                   Vec<2> (scales[patch])));
    return *(new (alloc) ScalarMappedElement<3> (local_ndof, order, basismat,
                                                 et, centers3[patch],
                                                 Vec<3> (scales[patch])));
  }

  template <int DIM>
  void AggregatedL2Space::UpdateGeometry (size_t patch, FlatArray<size_t> els)
  {
    Vec<DIM> center = 0.0;
    size_t npoints = 0;
    for (size_t el : els)
      for (int vertex : ma->GetElVertices (ElementId (VOL, el)))
        {
          center += ma->GetPoint<DIM> (vertex);
          npoints++;
        }
    center *= 1.0 / npoints;

    double diameter = 0.0;
    for (size_t el : els)
      for (int vertex : ma->GetElVertices (ElementId (VOL, el)))
        diameter
            = max (diameter, L2Norm (ma->GetPoint<DIM> (vertex) - center));

    if constexpr (DIM == 2)
      centers2[patch] = center;
    else
      centers3[patch] = center;
    scales[patch] = 1.0 / diameter;
  }

  template void
  AggregatedL2Space::UpdateGeometry<2> (size_t, FlatArray<size_t>);
  template void
  AggregatedL2Space::UpdateGeometry<3> (size_t, FlatArray<size_t>);
}

#ifdef NGS_PYTHON
void ExportAggregatedL2Space (py::module m)
{
  using namespace ngcomp;
  py::class_<AggregatedL2Space, shared_ptr<AggregatedL2Space>, FESpace> (
      m, "AggregatedL2Space")
      .def (py::init ([] (shared_ptr<MeshAccess> ma,
                          shared_ptr<ElementAggregation> aggregation,
                          py::kwargs kwargs) {
              Flags flags = ngcore::CreateFlagsFromKwArgs (kwargs);
              auto fes
                  = make_shared<AggregatedL2Space> (ma, flags, aggregation);
              fes->Update ();
              fes->FinalizeUpdate ();
              fes->ConnectAutoUpdate ();
              return fes;
            }),
            py::arg ("mesh"), py::arg ("aggregation"));
}
#endif
