#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>

#include "TrefftzElement.hpp"
#include "TrefftzFESpace.hpp"
#include "DiffOpMapped.hpp"

namespace ngcomp
{

  TrefftzFESpace :: TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "======== Constructor of TrefftzFESpace =========" << endl;
    cout << "Flags = " << flags << endl;

    order = int(flags.GetNumFlag ("order", 2));//flags.GetDefineFlag ("order");
		local_ndof = (BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1));

		int nel = ma->GetNE();
		ndof = local_ndof * nel;


    // needed for symbolic integrators and to draw solution
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<3,3>>>();
    //flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpMappedBoundary<3,3>>>();

    // (still) needed to draw solution
    //integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(),
    //                                             make_shared<ConstantCoefficientFunction>(1));
  }



  void TrefftzFESpace :: Update(LocalHeap & lh)
  {
		//int n_cell = ma->GetNE();
		//ndof = (BinCoeff(D + order, order) + BinCoeff(D + order-1, order-1)) * n_cell;
		cout << "update: order = " << order<< " ndof = " <<  ndof << endl;
	}

  void TrefftzFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
		int n_vert = ma->GetNV();
		int n_edge = ma->GetNEdges();
		int n_cell = ma->GetNE();
    // returns dofs of element ei
    // may be a volume triangle or boundary segment

		// returns dofs of element number elnr
		dnums.SetSize(0);

		Ngs_Element ngel = ma->GetElement (ei);

		// vertex dofs
		/*
		int local_nvert = 0;
		for (auto v : ngel.Vertices())
		{
			dnums.Append(v);
			local_nvert++;
		}

		for (int j = ei.Nr()*local_ndof+local_nvert; j-(ei.Nr()*local_ndof)<local_ndof; j++)
			{
					dnums.Append (j);
			}
*/
///////// TODO replace this by method from l2hofemspace
		for (int j = ei.Nr()*local_ndof; j-(ei.Nr()*local_ndof)<local_ndof; j++)
		{
			dnums.Append (j);
		}
	//cout << dnums;

	}


  FiniteElement & TrefftzFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
		return * new (alloc) TrefftzElement<3,3>;
  }

  /*
    register fe-spaces
    Object of type TrefftzFESpace can be defined in the pde-file via
    "define fespace v -type=trefftzfespace"
  */
  static RegisterFESpace<TrefftzFESpace> initi_trefftz ("trefftzfespace");
}



#ifdef NGS_PYTHON

void ExportTrefftzFESpace(py::module m)
{
  using namespace ngcomp;
  /*
    We just export the class here and use the FESpace constructor to create our space.
    This has the advantage, that we do not need to specify all the flags to parse (like
    dirichlet, definedon,...), but we can still append new functions only for that space.
   */
  py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>, FESpace>
    (m, "TrefftzFESpace", "FESpace with first order and second order trigs on 2d mesh")
    .def("Getndof", &TrefftzFESpace::Getndof);
}

#endif // NGS_PYTHON
