#ifndef FILE_AGGREGATEDL2SPACE_HPP
#define FILE_AGGREGATEDL2SPACE_HPP

#include "monomialfespace.hpp"
#include "dofpatches.hpp"
#include <xfem/aggregates.hpp>

namespace ngcomp
{
  class AggregatedL2Space : public MonomialFESpace, public DofPatchProvider
  {
    shared_ptr<ElementAggregation> aggregation;
    Array<int> element_to_patch;
    Array<Vec<2>> centers2;
    Array<Vec<3>> centers3;
    Array<double> scales;

  public:
    AggregatedL2Space (shared_ptr<MeshAccess> ama, const Flags &flags,
                       shared_ptr<ElementAggregation> aaggregation);

    void Update () override;
    void UpdateCouplingDofArray () override;
    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;
    DofPatches GetDofPatches () const override;

    string GetClassName () const override { return "AggregatedL2Space"; }
    shared_ptr<ElementAggregation> GetAggregation () const
    {
      return aggregation;
    }

  private:
    template <int DIM>
    void UpdateGeometry (size_t patch, FlatArray<size_t> els);
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportAggregatedL2Space (py::module m);
#endif

#endif
