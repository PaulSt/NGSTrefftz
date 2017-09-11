#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>

#include "TrefftzElement.hpp"
#include "TrefftzFESpace.hpp"

namespace ngcomp
{

  TrefftzFESpace ::TrefftzFESpace (shared_ptr<MeshAccess> ama,
                                   const Flags &flags)
      : FESpace (ama, flags)
  {
    this->D = 2;
    // int nel = ma->GetNE();
    // ndof = (BinCoeff(D + order, order) + BinCoeff(D + order-1, order-1)) *
    // nel;
    cout << "======== Constructor of TrefftzFESpace =========" << endl;
    cout << "Flags = " << flags << endl;

    order = flags.GetDefineFlag ("order");

    // needed for symbolic integrators and to draw solution
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>> ();
    flux_evaluator[VOL]
        = make_shared<T_DifferentialOperator<DiffOpGradient<2>>> ();
    evaluator[BND]
        = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>> ();

    // (still) needed to draw solution
    integrator[VOL] = GetIntegrators ().CreateBFI (
        "mass", ma->GetDimension (),
        make_shared<ConstantCoefficientFunction> (1));
  }

  void TrefftzFESpace ::Update (LocalHeap &lh)
  {

    // some global update:

    int n_vert = ma->GetNV ();
    int n_edge = ma->GetNEdges ();
    int n_cell = ma->GetNE ();

    first_edge_dof.SetSize (n_edge + 1);

    int ii = n_vert;
    for (int i = 0; i < n_edge; i++, ii += order - 1)
      first_edge_dof[i] = ii;
    first_edge_dof[n_edge] = ii;

    first_cell_dof.SetSize (n_cell + 1);
    for (int i = 0; i < n_cell; i++, ii += (order - 1) * (order - 2) / 2)
      first_cell_dof[i] = ii;
    first_cell_dof[n_cell] = ii;

    // cout << "first_edge_dof = " << endl << first_edge_dof << endl;
    // cout << "first_cell_dof = " << endl << first_cell_dof << endl;
    /*
                    ndof = ii;
    */
    ndof = (BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1))
           * n_cell;
  }

  void TrefftzFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment
    int locndof
        = BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1);

    // returns dofs of element number elnr
    dnums.SetSize (0);

    Ngs_Element ngel = ma->GetElement (ei);

    // vertex dofs
    for (auto v : ngel.Vertices ())
      dnums.Append (v);

    // edge dofs
    for (auto e : ngel.Edges ())
      {
        int first = first_edge_dof[e];
        int next = first_edge_dof[e + 1];
        for (int j = first; j < next; j++)
          dnums.Append (j);
      }
    if (ei.IsVolume ())
      {
        int first = first_cell_dof[ei.Nr ()];
        int next = first_cell_dof[ei.Nr () + 1];
        for (int j = first; j < next; j++)
          dnums.Append (j);
      }
  }

  FiniteElement &TrefftzFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    return *new (alloc) TrefftzElement<2, 3>;
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
  /*
    We just export the class here and use the FESpace constructor to create our
    space. This has the advantage, that we do not need to specify all the flags
    to parse (like dirichlet, definedon,...), but we can still append new
    functions only for that space.
   */
  py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>, FESpace> (
      m, "TrefftzFESpace",
      "FESpace with first order and second order trigs on 2d mesh")
      .def ("GetNVert", &TrefftzFESpace::GetNVert);
}

#endif // NGS_PYTHON
