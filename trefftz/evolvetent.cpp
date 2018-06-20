#include "evolvetent.hpp"
#include "trefftzfespace.hpp"
#include "tents/tents.hpp"
#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>

namespace ngcomp
{
  void EvolveTents (shared_ptr<MeshAccess> ma)
  {
    py::module et = py::module::import ("SolveTentSlab");
    auto tpmesh = NgsTPmesh (ma, 1);
    // Flags flag();
    // flag.SetFlag("order",4);
    // TrefftzFESpace fes(ma, Flags());
    // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
    // py::object ffes = py::cast(fes);
    // auto pyspace = py::class_<TrefftzFESpace,
    // shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str()); py::object
    // pyfes = et.attr("GetFESTrefftz")(ma); FESpace *ffes = pyfes.cast<FESpace
    // *>(); et.attr("EvolveTent")(pyfes,?,?);
  }

  shared_ptr<MeshAccess> NgsTPmesh (shared_ptr<MeshAccess> ma, float wavespeed)
  {
    Point2IndexMap *pim = new Point2IndexMap (); // Your map type may vary,
                                                 // just change the typedef

    int basedim = ma->GetDimension ();
    TentPitchedSlab<1> tps
        = TentPitchedSlab<1> (ma); // collection of tents in timeslab
    if (basedim >= 2)
      cout << "oh no" << endl;
    tps.PitchTents (1.0, wavespeed); // adt = time slab height, wavespeed

    auto mesh = make_shared<netgen::Mesh> ();
    mesh->SetDimension (ma->GetDimension () + 1);
    netgen::SetGlobalMesh (mesh); // for visualization
    mesh->SetGeometry (make_shared<netgen::NetgenGeometry> ());

    Vec<1> point;
    netgen::FaceDescriptor *fd = new netgen::FaceDescriptor ();
    fd->SetSurfNr (1);
    fd->SetDomainIn (1);
    fd->SetDomainOut (0);
    fd->SetBCProperty (1);
    mesh->AddFaceDescriptor (*fd);
    for (Tent *tent : tps.tents)
      {
        Vector<netgen::PointIndex> vertices (tent->nbv.Size () + 2);
        point = ma->GetPoint<1> (tent->vertex);
        vertices[tent->nbv[0] > tent->vertex ? 0 : 2] = AddPointUnique (
            mesh, pim, netgen::Point3d (point (0), tent->tbot, 0));
        vertices[tent->nbv[0] > tent->vertex ? 2 : 0] = AddPointUnique (
            mesh, pim, netgen::Point3d (point (0), tent->ttop, 0));
        for (int k = 0; k < tent->nbv.Size (); k++)
          {
            point = ma->GetPoint<1> (tent->nbv[k]);
            vertices[2 * k + 1] = AddPointUnique (
                mesh, pim, netgen::Point3d (point (0), tent->nbtime[k], 0));
          }
        int index = 1;
        netgen::Element2d *newel = nullptr;
        if (tent->nbv.Size () == 1)
          {
            newel = new netgen::Element2d (netgen::TRIG);
            for (int i = 0; i < 3; i++)
              (*newel)[i] = vertices[i];
            newel->SetIndex (index);
          }
        else if (tent->nbv.Size () == 2)
          {
            newel = new netgen::Element2d (netgen::QUAD);
            for (int i = 0; i < 4; i++)
              (*newel)[i] = vertices[i];
            newel->SetIndex (index);
          }
        mesh->AddSurfaceElement (*newel);
        // cout << *tent << endl;
      }
    auto tpmesh = make_shared<MeshAccess> (mesh);
    return tpmesh;
  }

  netgen::PointIndex Point2Index (Point2IndexMap *pim, netgen::Point3d p)
  {
    Point2IndexMap::iterator lb = pim->find (p);

    if (lb != pim->end () && !(pim->key_comp () (p, lb->first)))
      return lb->second;
    else
      {
        // the key does not exist in the map
        // add it to the map
        netgen::PointIndex newpind (pim->size () + 1);
        pim->insert (lb, Point2IndexMap::value_type (
                             p, newpind)); // Use lb as a hint to insert,
        return newpind;
      }
  }

  netgen::PointIndex AddPointUnique (shared_ptr<netgen::Mesh> ngma,
                                     Point2IndexMap *pim, netgen::Point3d p)
  {
    netgen::PointIndex pi = Point2Index (pim, p);
    if (pi == pim->size ())
      ngma->AddPoint (p);
    return pi;
  }

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m)
{
  m.def (
      "EvolveTents", [] (shared_ptr<MeshAccess> ma) //-> shared_ptr<MeshAccess>
      { return EvolveTents (ma); } //, py::call_guard<py::gil_scoped_release>()
  );
  m.def ("NgsTPmesh",
         [] (shared_ptr<MeshAccess> ma,
             float wavespeed) -> shared_ptr<MeshAccess> {
           return NgsTPmesh (ma, wavespeed);
         } //, py::call_guard<py::gil_scoped_release>()
  );
}
#endif // NGS_PYTHON
