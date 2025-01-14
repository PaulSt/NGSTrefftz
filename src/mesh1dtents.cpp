#include <tents.hpp>
#include "mesh1dtents.hpp"

namespace ngcomp
{
  typedef map<netgen::Point3d, netgen::PointIndex> Point2IndexMap;

  netgen::PointIndex Point2Index (Point2IndexMap *pim, netgen::Point3d p)
  {
    Point2IndexMap::iterator lb = pim->find (p);

    if (lb != pim->end () && !(pim->key_comp () (p, lb->first)))
      return lb->second;
    else
      {
        // the key does not exist in the map add it to the map
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
    if (size_t (pi) == pim->size ())
      ngma->AddPoint (p);
    return pi;
  }

  shared_ptr<MeshAccess> NgsTP1dmesh (shared_ptr<TentPitchedSlab> tps)
  {
    Point2IndexMap *pim = new Point2IndexMap (); // Your map type may vary,
                                                 // just change the typedef
    int index = 1;
    double dt_eps = 0.00001;
    shared_ptr<MeshAccess> ma = tps->ma;
    double dt = tps->GetSlabHeight ();
    // get boundaries of the init mesh
    Array<double> bd_points (0);
    for (auto el : ma->Elements (BND))
      bd_points.Append (ma->GetPoint<1> (el.Vertices ()[0]) (0));

    auto mesh = make_shared<netgen::Mesh> ();
    mesh->SetDimension (ma->GetDimension () + 1);
    netgen::SetGlobalMesh (mesh); // for visualization
    mesh->SetGeometry (make_shared<netgen::NetgenGeometry> ());
    mesh->SetMaterial (1, "mat");

    netgen::FaceDescriptor fd (1, 1, 0, 1);
    // int ind_fd =
    mesh->AddFaceDescriptor (fd);

    netgen::FaceDescriptor fdi (1, 1, 0, 1);
    fdi.SetBCName (new string ("inflow"));
    int ind_fdi = mesh->AddFaceDescriptor (fdi);

    netgen::FaceDescriptor fdo (2, 1, 0, 2);
    fdo.SetBCName (new string ("outflow"));
    int ind_fdo = mesh->AddFaceDescriptor (fdo);

    netgen::FaceDescriptor fdd (3, 1, 0, 3);
    fdd.SetBCName (new string ("dirichlet"));
    int ind_fdd = mesh->AddFaceDescriptor (fdd);

    for (int i = 0; i < tps->GetNTents (); i++)
      {
        const Tent *tent = &tps->GetTent (i);
        // cout << *tent << endl;
        if (tent->ttop < tent->tbot + dt_eps)
          {
            continue;
            cout << "had to skip degenerate tent" << endl;
          }
        // Add vertices and 2d Elements to the mesh
        Vector<netgen::PointIndex> vertices (tent->nbv.Size () + 2);
        double pointc = ma->GetPoint<1> (tent->vertex) (0);
        int ibot = tent->nbv[0] > tent->vertex ? 0 : 2;
        int itop = tent->nbv[0] > tent->vertex ? 2 : 0;
        vertices[ibot] = AddPointUnique (
            mesh, pim, netgen::Point3d (pointc, tent->tbot, 0));
        vertices[itop] = AddPointUnique (
            mesh, pim, netgen::Point3d (pointc, tent->ttop, 0));
        for (size_t k = 0; k < tent->nbv.Size (); k++)
          {
            vertices[2 * k + 1] = AddPointUnique (
                mesh, pim,
                netgen::Point3d (ma->GetPoint<1> (tent->nbv[k]) (0),
                                 tent->nbtime[k], 0));
          }
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

        for (size_t k = 0; k < tent->nbv.Size (); k++)
          {
            // Add 1d Elements - inflow
            if (tent->tbot < dt_eps && tent->nbtime[k] < dt_eps)
              {
                netgen::Segment *newel = new netgen::Segment ();
                (*newel)[tent->vertex < tent->nbv[k] ? 0 : 1] = vertices[ibot];
                (*newel)[tent->vertex < tent->nbv[k] ? 1 : 0]
                    = vertices[2 * k + 1];
                newel->si = ind_fdi;
                newel->edgenr = index;
                newel->epgeominfo[0].edgenr = index;
                newel->epgeominfo[1].edgenr = index;
                mesh->AddSegment (*newel);
              }

            // Add 1d Elements - outflow
            if (tent->ttop > dt - dt_eps && tent->nbtime[k] > dt - dt_eps)
              {
                netgen::Segment *newel = new netgen::Segment ();
                (*newel)[tent->vertex < tent->nbv[k] ? 1 : 0] = vertices[itop];
                (*newel)[tent->vertex < tent->nbv[k] ? 0 : 1]
                    = vertices[2 * k + 1];
                newel->si = ind_fdo;
                newel->edgenr = index;
                newel->epgeominfo[0].edgenr = index;
                newel->epgeominfo[1].edgenr = index;
                mesh->AddSegment (*newel);
              }

            // Add 1d Elements - dirichlet
            if (pointc == bd_points[0] || pointc == bd_points[1])
              {
                netgen::Segment *newel = new netgen::Segment ();
                (*newel)[0] = vertices[2];
                (*newel)[1] = vertices[0];
                newel->si = ind_fdd;
                newel->edgenr = index;
                newel->epgeominfo[0].edgenr = index;
                newel->epgeominfo[1].edgenr = index;
                mesh->AddSegment (*newel);
              }
          }
      }

    mesh->SetBCName (1, "inflow");
    mesh->SetBCName (2, "outflow");
    mesh->SetBCName (3, "dirichlet");

    auto tpmesh = make_shared<MeshAccess> (mesh);
    // cout << tpmesh->GetMaterials(BND) << endl;
    return tpmesh;
  }

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMesh1dTents (py::module m)
{
  m.def ("Mesh1dTents",
         [] (shared_ptr<TentPitchedSlab> tps) -> shared_ptr<MeshAccess> {
           return NgsTP1dmesh (tps);
         } //, py::call_guard<py::gil_scoped_release>()
  );
}
#endif // NGS_PYTHON
