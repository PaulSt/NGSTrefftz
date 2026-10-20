#include <comp.hpp> // provides FESpace, ...
#include <python_comp.hpp>

#include "pufe.hpp"
#include "pufespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{
  template <int D, typename OP>
  class DiffOpPUEnrichment : public DiffOp<DiffOpPUEnrichment<D, OP>>
  {
  public:
    enum
    {
      DIM = OP::DIM,
      DIM_SPACE = D,
      DIM_ELEMENT = D,
      DIM_DMAT = OP::DIM_DMAT,
      DIFFORDER = OP::DIFFORDER
    };
    static auto GetDimensions () { return OP::GetDimensions (); }
    static string Name () { return "enrichment_" + OP::Name (); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement &fel, const MIP &mip,
                                MAT &&mat, LocalHeap &lh)
    {
      auto &pu = static_cast<const PUFElement<D> &> (fel);
      int n = fel.GetNDof () / (D + 1);
      for (int vertex = 0; vertex < D + 1; vertex++)
        OP::GenerateMatrix (pu.GetEnrichment (vertex), mip,
                            mat.Cols (vertex * n, (vertex + 1) * n), lh);
    }
  };

  PUFESpace::PUFESpace (shared_ptr<MeshAccess> ama, const Flags &flags, bool)
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
          additional_evaluators.Set (
              "enrichment", make_shared<T_DifferentialOperator<
                                DiffOpPUEnrichment<2, DiffOpMapped<2>>>> ());
          additional_evaluators.Set (
              "enrichment_hesse",
              make_shared<T_DifferentialOperator<
                  DiffOpPUEnrichment<2, DiffOpMappedHesse<2>>>> ());
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

  void PUFESpace::Update ()
  {
    this->local_ndof = BinCoeff (D + order, order);
    int nvert = ma->GetNV ();

    size_t ndof = local_ndof * nvert;
    SetNDof (ndof);

    // FESpace::Update ();
    UpdateCouplingDofArray ();
  }

  void PUFESpace::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
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

  DofPatches PUFESpace::GetDofPatches () const
  {
    TableCreator<size_t> elements (ma->GetNV ());
    for (; !elements.Done (); elements++)
      for (size_t vertex : Range (ma->GetNV ()))
        for (size_t el : ma->GetVertexElements (vertex))
          if (DefinedOn (ElementId (VOL, el)))
            elements.Add (vertex, el);
    DofPatches patches{ elements.MoveTable (),
                        Table<DofId> (ma->GetNV (), local_ndof) };
    auto dofs = patches.dofs.AsArray ();
    for (size_t d : Range (dofs))
      dofs[d] = d;
    return patches;
  }

  void PUFESpace::UpdateCouplingDofArray ()
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

  FiniteElement &PUFESpace::GetFE (ElementId ei, Allocator &alloc) const
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

  DocInfo PUFESpace::GetDocu ()
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
