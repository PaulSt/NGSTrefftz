#include "evolvetent.hpp"
#include "trefftzfespace.hpp"
#include "tents/tents.hpp"
#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>


namespace ngcomp
{
  void EvolveTents(shared_ptr<MeshAccess> ma)
  {
    py::module et = py::module::import("SolveTentSlab");
    auto tpmesh = ngs_tpmesh(ma, 1);
    // Flags flag();
    // flag.SetFlag("order",4);
    // TrefftzFESpace fes(ma, Flags());
    // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
    // py::object ffes = py::cast(fes);
    // auto pyspace = py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str());
    // py::object pyfes = et.attr("GetFESTrefftz")(ma);
    // FESpace *ffes = pyfes.cast<FESpace *>();
    // et.attr("EvolveTent")(pyfes,?,?);
  }


  shared_ptr<MeshAccess> ngs_tpmesh(shared_ptr<MeshAccess> ma, float wavespeed)
  {
    map<netgen::Point3d, netgen::PointIndex> point2index_map;    // Your map type may vary, just change the typedef

    cout << point2index(&point2index_map, netgen::Point3d(1,1,1)) << endl;
    cout << point2index(&point2index_map, netgen::Point3d(2,1,1)) << endl;
    cout << point2index(&point2index_map, netgen::Point3d(1,1,1)) << endl;
    cout << point2index(&point2index_map, netgen::Point3d(1,2,1)) << endl;
    cout << point2index(&point2index_map, netgen::Point3d(2,1,1)) << endl;

    int basedim =  ma->GetDimension(); 
    TentPitchedSlab<1> tps = TentPitchedSlab<1>(ma);      // collection of tents in timeslab
    if(basedim >= 2)
      cout << "oh no" << endl;
    tps.PitchTents(1.0, wavespeed); // adt = time slab height, wavespeed

    auto mesh = make_shared<netgen::Mesh>();
    mesh -> SetDimension(ma->GetDimension() + 1);
    netgen::SetGlobalMesh(mesh);  // for visualization
    mesh -> SetGeometry (make_shared<netgen::NetgenGeometry>());

    Vec<1> point;
    netgen::FaceDescriptor * fd = new netgen::FaceDescriptor();
    fd->SetSurfNr(1);
    fd->SetDomainIn(1);
    fd->SetDomainOut(0);
    fd->SetBCProperty(1);
    mesh -> AddFaceDescriptor (*fd);
    for(Tent* tent : tps.tents)
    { 
      Vector<netgen::PointIndex> vertices(tent->nbv.Size()+2);
      point = ma -> GetPoint<1>(tent->vertex); 
      vertices[ tent->nbv[0]>tent->vertex ? 0 : 2 ] = mesh -> AddPoint(netgen::Point3d(point(0),tent->tbot,0));
      vertices[ tent->nbv[0]>tent->vertex ? 2 : 0 ] = mesh -> AddPoint(netgen::Point3d(point(0),tent->ttop,0));
      for (int k = 0; k < tent->nbv.Size(); k++)
      {
        point = ma->GetPoint<1>(tent->nbv[k]);
        vertices[2*k+1] = mesh -> AddPoint(netgen::Point3d(point(0),tent->nbtime[k],0));
      }
      int index = 1;
      netgen::Element2d * newel = nullptr;
      if (tent->nbv.Size() == 1)
      {
       newel = new netgen::Element2d(netgen::TRIG);
       for (int i = 0; i < 3; i++)
         (*newel)[i] = vertices[i];
       newel->SetIndex(index);
      }
      else if (tent->nbv.Size() == 2)
      {
       newel = new netgen::Element2d(netgen::QUAD);
       for (int i = 0; i < 4; i++)
         (*newel)[i] = vertices[i];
       newel->SetIndex(index);
      }
      mesh -> AddSurfaceElement(*newel);
      // cout << *tent << endl;
    }
    auto tpmesh = make_shared<MeshAccess>(mesh);
    return tpmesh;
  }

  netgen::PointIndex point2index(map<netgen::Point3d, netgen::PointIndex> *point2index_map, netgen::Point3d p)
  {
    map<netgen::Point3d, netgen::PointIndex>::iterator lb = point2index_map->find(p);

    if(lb != point2index_map->end() && !(point2index_map->key_comp()(p, lb->first)))
    {
      cout << "found point " << lb->first << " with pind: " << lb->second << endl;
      return lb->second;
    }
    else
    {
      // the key does not exist in the map
      // add it to the map
      netgen::PointIndex newpind(0);
      if(lb == point2index_map->begin())
      {
        cout << "adding first point index" << endl;
      }
      else
      {
        lb--;
        newpind = lb->second;
        lb++;
      }
      cout << "insterting: " << p << " with " << newpind+1 << endl;
      point2index_map->insert(lb, map<netgen::Point3d, netgen::PointIndex>::value_type(p, ++newpind));    // Use lb as a hint to insert,
      return newpind;
    }
  }
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m)
{
	m.def("EvolveTents", [](shared_ptr<MeshAccess> ma) //-> shared_ptr<MeshAccess>
					{
            return EvolveTents(ma);
					}//, py::call_guard<py::gil_scoped_release>()
      );
	m.def("ngs_tpmesh", [](shared_ptr<MeshAccess> ma, float wavespeed) -> shared_ptr<MeshAccess>
					{
            return ngs_tpmesh(ma,wavespeed);
					}//, py::call_guard<py::gil_scoped_release>()
      );
}
#endif // NGS_PYTHON

