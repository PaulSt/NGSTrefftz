#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>

#include "TrefftzElement.hpp"
#include "TrefftzFESpace.hpp"
#include "DiffOpMapped.hpp"


namespace ngcomp
{

  TrefftzFESpace :: TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
		DefineNumFlag("wavespeed");
    cout << "======== Constructor of TrefftzFESpace =========" << endl;
    cout << "Flags = " << flags;

		D = ma->GetDimension();

    order = int(flags.GetNumFlag ("order", 3));//flags.GetDefineFlag ("order");
		c = flags.GetNumFlag ("wavespeed", 1);

		local_ndof = (BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1));
		int nel = ma->GetNE();
		ndof = local_ndof * nel;

		switch (D) {
			case 2:
			{
		    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<2,TrefftzElement<2,3>>>>();
				flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<2, TrefftzElement<2,3>>>>();
				evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpMappedBoundary<2,TrefftzElement<2,3>>>>();
				break;
			}
			case 3:
			{
				// needed for symbolic integrators and to draw solution
		    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<3,TrefftzElement<3,3>>>>();
		    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<3, TrefftzElement<3,3>>>>();
		    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpMappedBoundary<3, TrefftzElement<3,3>>>>();
		    // (still) needed to draw solution
		    //integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(),
		    //                                             make_shared<ConstantCoefficientFunction>(1));
				break;
			}
		}
  }



  void TrefftzFESpace :: Update(LocalHeap & lh)
  {
		local_ndof = (BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1));
		int nel = ma->GetNE();
		ndof = local_ndof * nel;

		cout << "update: order = " << order << " D: " << D << " ndof = " <<  ndof << " local_ndof:" << local_ndof << endl <<
		"================================================" << endl ;
	}

  void TrefftzFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
		int n_vert = ma->GetNV();
		int n_edge = ma->GetNEdges();
		int n_cell = ma->GetNE();
		dnums.SetSize(0);
		Ngs_Element ngel = ma->GetElement (ei);
		for (int j = ei.Nr()*local_ndof; j-(ei.Nr()*local_ndof)<local_ndof; j++)
		{
			dnums.Append (j);
		}
	//cout << dnums;
	}


  FiniteElement & TrefftzFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
		auto vertices_index = ma->GetElVertices(ei);
		//cout << "element vectice coord: \n"  << ma->GetPoint<3>(vertices_index[0]) << endl<< ma->GetPoint<3>(vertices_index[1]) <<endl<<ma->GetPoint<3>(vertices_index[2])<<endl<<ma->GetPoint<3>(vertices_index[3])<<endl;
		if(order != 3){cout << "order not yet supported"<<endl;}
		switch (D) {
			case 2:
			{
				Vec<2> center = 0;
				for(int d=0;d<3;d++) center += ma->GetPoint<2>(vertices_index[d]);
				center *= (1.0/3.0);
				return *(new (alloc) TrefftzElement<2,3>) ->SetCenter(center) ->SetWavespeed(c) ->SetElSize(1);
				break;
			}
			case 3:
			{
				Vec<3> center = 0;
				for(int d=0;d<4;d++) center += ma->GetPoint<3>(vertices_index[d]);
				center *= 0.25;
				return *(new (alloc) TrefftzElement<3,3>) ->SetCenter(center) ->SetWavespeed(c) ->SetElSize(1);
				break;
			}
		}
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
	using namespace ngfem;
  /*
    We just export the class here and use the FESpace constructor to create our space.
    This has the advantage, that we do not need to specify all the flags to parse (like
    dirichlet, definedon,...), but we can still append new functions only for that space.
   */
  py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>, FESpace>
    (m, "TrefftzFESpace", "FESpace with first order and second order trigs on 2d mesh")
    .def("GetNDof", &TrefftzFESpace::GetNDof)
		;
}

#endif // NGS_PYTHON
