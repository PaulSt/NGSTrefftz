#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>

#include "trefftzwavefe.hpp"
#include "trefftzfespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{

  TrefftzFESpace ::TrefftzFESpace (shared_ptr<MeshAccess> ama,
                                   const Flags &flags)
      : FESpace (ama, flags)
  {

    DefineNumFlag ("wavespeed");
    cout << "======== Constructor of TrefftzFESpace =========" << endl;
    cout << "Flags:" << endl << flags;

    D = ma->GetDimension ();

    order = int (
        flags.GetNumFlag ("order", 3)); // flags.GetDefineFlag ("order");
    c = flags.GetNumFlag ("wavespeed", 1);
    basistype = flags.GetNumFlag ("basistype", 0);
    useshift = flags.GetNumFlag ("useshift", 1);

    local_ndof = (BinCoeff (D - 1 + order, order)
                  + BinCoeff (D - 1 + order - 1, order - 1));
    nel = ma->GetNE ();
    ndof = local_ndof * nel;

    switch (D)
      {
      case 1:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<1>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<1>>> ();
          break;
        }
      case 2:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<2>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<2>>> ();
          break;
        }
      case 3:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<3>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<3>>> ();
          break;
        }
      }
  }

  void TrefftzFESpace ::Update (LocalHeap &lh)
  {
    FESpace::Update (lh);
    cout << "update: order = " << order << " D: " << D << " ndof = " << ndof
         << " local_ndof:" << local_ndof << endl
         << "================================================" << endl;
  }

  void TrefftzFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize (0);
    if (!DefinedOn (ei) || ei.VB () != VOL)
      return;
    // int n_vert = ma->GetNV();		int n_edge = ma->GetNEdges();		int n_cell =
    // ma->GetNE(); Ngs_Element ngel = ma->GetElement (ei);
    for (int j = ei.Nr () * local_ndof; j < local_ndof * (ei.Nr () + 1); j++)
      {
        dnums.Append (j);
      }
    // cout << "GetDofNrs: ei.Nr() = " << ei.Nr() << " local_ndof:" <<
    // local_ndof << " ndof: " << ndof << " dnums: \n" << dnums << endl <<
    // "================================================" << endl ;
  }

  FiniteElement &TrefftzFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    auto vertices_index = ma->GetElVertices (ei);
    // cout << "element vectice coord: \n"  <<
    // ma->GetPoint<3>(vertices_index[0]) << endl<<
    // ma->GetPoint<3>(vertices_index[1])
    // <<endl<<ma->GetPoint<3>(vertices_index[2])<<endl<<ma->GetPoint<3>(vertices_index[3])<<endl;
    switch (ma->GetElType (ei))
      {
      case ET_SEGM:
        {
          return *(new (alloc) TrefftzWaveFE<1> (order, c, ET_SEGM, basistype))
                      ->SetCenter (ElCenter<1> (ei))
                      ->SetElSize (Adiam<1> (ei));
          break;
        }
      case ET_QUAD:
      case ET_TRIG:
        {
          return *(new (alloc) TrefftzWaveFE<2> (order, c, ma->GetElType (ei),
                                                 basistype))
                      ->SetCenter (ElCenter<2> (ei))
                      ->SetElSize (Adiam<2> (ei));
          break;
        }
      case ET_HEX:
      case ET_PRISM:
      case ET_PYRAMID:
      case ET_TET:
        {
          return *(new (alloc) TrefftzWaveFE<3> (order, c, ma->GetElType (ei),
                                                 basistype))
                      ->SetCenter (ElCenter<3> (ei))
                      ->SetElSize (Adiam<3> (ei));
          break;
        }
      }
    return *(new (alloc) TrefftzWaveFE<1> ());
  }

  template <int D> double TrefftzFESpace ::Adiam (ElementId ei) const
  {
    double anisotropicdiam = 0.0;
    auto vertices_index = ma->GetElVertices (ei);
    for (auto vertex1 : vertices_index)
      {
        for (auto vertex2 : vertices_index)
          {
            Vec<D> v1 = ma->GetPoint<D> (vertex1);
            Vec<D> v2 = ma->GetPoint<D> (vertex2);
            // cout << "v1: " << v1 << " v1 part: " << v1(1,D-1) << "norm " <<
            // L2Norm(v1) << endl ;
            anisotropicdiam
                = max (anisotropicdiam,
                       sqrt (L2Norm2 (v1 (0, D - 2) - v2 (0, D - 2))
                             + pow (c * (v1 (D - 1) - v2 (D - 1)), 2)));
          }
      }
    return anisotropicdiam * useshift + 1 * (useshift == 0);
  }

  template <int D> Vec<D> TrefftzFESpace ::ElCenter (ElementId ei) const
  {
    Vec<D> center = 0;
    auto vertices_index = ma->GetElVertices (ei);
    for (auto vertex : vertices_index)
      center += ma->GetPoint<D> (vertex);
    center *= (1.0 / vertices_index.Size ()) * useshift;
    return center;
  }

  /*
     register fe-spaces
     Object of type TrefftzFESpace can be defined in the pde-file via
     "define fespace v -type=trefftzfespace"
     */
  static RegisterFESpace<TrefftzFESpace> initi_trefftz ("trefftzfespace");
}

#ifdef NGS_PYTHON

void ExportTrefftzFESpace (py::module m)
{
  using namespace ngcomp;
  using namespace ngfem;
  /*
     We just export the class here and use the FESpace constructor to create
     our space. This has the advantage, that we do not need to specify all the
     flags to parse (like dirichlet, definedon,...), but we can still append
     new functions only for that space.
     */
  py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>, FESpace> (
      m, "TrefftzFESpace",
      "FESpace with first order and second order trigs on 2d mesh")
      .def ("GetNDof", &TrefftzFESpace::GetNDof);
  m.def ("GetNDof", [] (shared_ptr<FESpace> fes) {
    cout << typeid (*fes).name () << endl;
    // fes->GetNDof();
  });
}

#endif // NGS_PYTHON
