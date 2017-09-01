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
    cout << "Update TrefftzFESpace, #vert = " << ma->GetNV ()
         << ", #edge = " << ma->GetNEdges ()
         << ", dim: " << ma->GetDimension () << endl;

    this->D = 2;
    int nel = ma->GetNE ();
    cout << "================= NEL" << nel << endl;
    ndof = (BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1))
           * nel;
    cout << "#dof: " << ndof << endl;
  }

  void TrefftzFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment

    int locndof
        = BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1);
    // cout << locndof;
    dnums.SetSize (0);
    /*
        // first dofs are vertex numbers:
        for (auto v : ma->GetElVertices(ei))
          dnums.Append (v);

            for (auto e : ma->GetElEdges(ei))
          dnums.Append (nvert+e);
                            */
  }

  FiniteElement &TrefftzFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    return *new (alloc) TrefftzElement<2, 4>;
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
