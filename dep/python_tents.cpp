#include "tents.hpp"
#include <python_ngstd.hpp>

typedef CoefficientFunction CF;

// python export of tent mesh
auto ExportTimeSlab(py::module &m)
{
  auto pyname = "TentSlab";
  auto pydocu = "Tent pitched slab in D + 1 time dimensions";
  py::class_<TentPitchedSlab, shared_ptr<TentPitchedSlab>>
    (m, pyname, pydocu)
    .def(py::init([](shared_ptr<MeshAccess> ma, string method_name, int heapsize)
		  {
		    ngstents::PitchingMethod method = [method_name]{
		      if(method_name == "edge") return ngstents::EEdgeGrad;
		      else if(method_name == "vol") return ngstents::EVolGrad;
		      else//just for static analyzers. the code should not reach this case
			{
			  cout << "Invalid method! Setting edge algorithm as default..." << endl;
			  return ngstents::EEdgeGrad;
			}
		    }();
		    auto tps = TentPitchedSlab(ma,  heapsize);
		    tps.SetPitchingMethod(method);
		    return tps;
		  }),
	 py::arg("mesh"), py::arg("method") = "edge", py::arg("heapsize") = 1000000
	 )
    .def_readonly("mesh", &TentPitchedSlab::ma)
    .def_property_readonly("gradphi",[](shared_ptr<TentPitchedSlab> self) -> shared_ptr<CF>
			   {
			     return self->cfgradphi;
			   })
    .def("SetMaxWavespeed", [](shared_ptr<TentPitchedSlab> self, py::object wavespeed)
	 {
	   if (auto ws = py::extract<double> (wavespeed); ws.check())
	     self->SetMaxWavespeed(ws());
           else if (auto ws = py::extract<shared_ptr<CF>>(wavespeed); ws.check())
             self->SetMaxWavespeed(ws());
           else
             throw Exception("wrong argument type in SetMaxWavespeed");
         })
    .def("PitchTents",[](shared_ptr<TentPitchedSlab> self,
			 const double dt, const bool local_ct, const double global_ct)
	 {
	   int dim = self->ma->GetDimension();
	   bool success = false;
	   switch(dim){
	   case 1:
	     success = self->PitchTents<1>(dt,local_ct,global_ct);
	     break;
	   case 2:
	     success = self->PitchTents<2>(dt,local_ct,global_ct);
	     break;
	   case 3:
	     success = self->PitchTents<3>(dt,local_ct,global_ct);
	     break;
	   default:
	     throw Exception("TentPitchedSlab not avaiable for dimension "+ToString(dim));
	   }
	   return success;
	 },
	 py::arg("dt"), py::arg("local_ct") = false, py::arg("global_ct") = 1.0)
    .def("GetNTents", &TentPitchedSlab::GetNTents)
    .def("GetNLayers", &TentPitchedSlab::GetNLayers)
    .def("GetSlabHeight", &TentPitchedSlab::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab::MaxSlope)
    .def("GetTent", &TentPitchedSlab::GetTent, pybind11::return_value_policy::reference_internal)
    ////////////////////////////
    // visualization functions
    ////////////////////////////
    .def("_TentData1D", [](shared_ptr<TentPitchedSlab> self)
	 {
	   py::list ret;
	   for(int i = 0; i < self->GetNTents(); i++)
	     {
	       const Tent & tent = self->GetTent(i);
	       py::list reti;
	       reti.append(py::make_tuple(tent.vertex, tent.ttop,
					  tent.tbot, tent.level));
	       for(int j = 0; j< tent.nbv.Size(); j++)
		 reti.append(py::make_tuple(tent.nbv[j],tent.nbtime[j]));
	       ret.append(reti);
	     }
	   return ret;
	 })
    .def("DrawPitchedTentsVTK", [](shared_ptr<TentPitchedSlab> self, string vtkfilename)
	 {
	   if(self->ma->GetDimension() != 2)
	     throw Exception("VTK export is only supported for 2D spatial meshes");
	   else
	     self->DrawPitchedTentsVTK(vtkfilename);
	 },
	 py::arg("vtkfilename") = "output")
    .def("DrawPitchedTentsGL", [](shared_ptr<TentPitchedSlab> self)
  	 {
	   if(self->ma->GetDimension() == 1)
	     throw Exception("Not supported for 1D spatial meshes");
	   
  	   int nlevels;
  	   Array<int> tentdata;
  	   Array<double> tenttimes;
  	   self->DrawPitchedTentsGL(tentdata, tenttimes, nlevels);
  	   py::list data, times;
  	   for(auto i : Range(tentdata))
  	     {
  	       data.append(tentdata[i]);
  	       // note: time values make sense only in 2D case.
  	       // They are not used in 3D case, i.e. they are
  	       // ignored by tents_visualization (ngsgui) and webgui.
  	       times.append(tenttimes[i]);
  	     }
  	   return py::make_tuple(data,times,self->GetNTents(),nlevels);
  	 })
    ;
}


void ExportTents(py::module & m) {

  py::class_<Tent, shared_ptr<Tent>>(m, "Tent", "Tent structure")
    .def_readonly("vertex", &Tent::vertex)
    .def_readonly("ttop", &Tent::ttop)
    .def_readonly("tbot", &Tent::tbot)
    .def_readonly("nbv", &Tent::nbv)
    .def_readonly("nbtime", &Tent::nbtime)
    .def_readonly("els", &Tent::els)
    .def_readonly("level", &Tent::level)
    .def_readonly("internal_facets", &Tent::internal_facets)
    .def("MaxSlope", &Tent::MaxSlope);

  ExportTimeSlab(m);
}


PYBIND11_MODULE(_pytents, m) {
  py::module_::import("ngsolve");
  m.attr("__name__") = "ngstents";
  m.attr("__package__") = "ngstents";
  ExportTents(m);
}
