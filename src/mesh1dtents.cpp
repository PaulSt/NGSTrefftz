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

    // the key does not exist in the map add it to the map
    netgen::PointIndex newpind (pim->size () + 1);
    pim->insert (lb, Point2IndexMap::value_type (p, newpind));
    return newpind;
  }

  netgen::PointIndex AddPointUnique (shared_ptr<netgen::Mesh> ngma,
                                     Point2IndexMap *pim, netgen::Point3d p)
  {
    const size_t old_size = pim->size ();
    netgen::PointIndex pi = Point2Index (pim, p);
    if (pim->size () != old_size)
      ngma->AddPoint (p);
    return pi;
  }

  int AddBoundaryEdgeDescriptor (netgen::Mesh &mesh, int edge_nr,
                                 int face_descriptor_index, const string &name)
  {
    netgen::EdgeDescriptor ed;
    ed.SetEdgeNr (edge_nr);
    ed.SetName (name);
    ed.SetIndex (face_descriptor_index);
    ed.SetSurfNr (0, 1);
    ed.SetSurfNr (1, 0);
    ed.SetDomainIn (1);
    ed.SetDomainOut (0);
    ed.SetTLOSurface (edge_nr);
    return mesh.AddEdgeDescriptor (ed);
  }

  shared_ptr<MeshAccess> NgsTP1dmesh (shared_ptr<TentPitchedSlab> tps)
  {
    Point2IndexMap pim;
    constexpr int domain_index = 1;
    constexpr double dt_eps = 0.00001;

    shared_ptr<MeshAccess> ma = tps->ma;
    const double dt = tps->GetSlabHeight ();

    Array<double> bd_points (0);
    for (auto el : ma->Elements (BND))
      bd_points.Append (ma->GetPoint<1> (el.Vertices ()[0]) (0));

    auto mesh = make_shared<netgen::Mesh> ();
    mesh->SetDimension (ma->GetDimension () + 1);
    netgen::SetGlobalMesh (mesh); // for visualization
    mesh->SetGeometry (make_shared<netgen::NetgenGeometry> ());
    mesh->SetMaterial (domain_index, "mat");

    netgen::FaceDescriptor fd (1, 1, 0, 1);
    mesh->AddFaceDescriptor (fd);

    netgen::FaceDescriptor fdi (1, 1, 0, 1);
    fdi.SetBCName ("inflow");
    const int ind_fdi = mesh->AddFaceDescriptor (fdi);

    netgen::FaceDescriptor fdo (2, 1, 0, 2);
    fdo.SetBCName ("outflow");
    const int ind_fdo = mesh->AddFaceDescriptor (fdo);

    netgen::FaceDescriptor fdd (3, 1, 0, 3);
    fdd.SetBCName ("dirichlet");
    const int ind_fdd = mesh->AddFaceDescriptor (fdd);

    const int ed_inflow
        = AddBoundaryEdgeDescriptor (*mesh, 1, ind_fdi, "inflow");
    const int ed_outflow
        = AddBoundaryEdgeDescriptor (*mesh, 2, ind_fdo, "outflow");
    const int ed_dirichlet
        = AddBoundaryEdgeDescriptor (*mesh, 3, ind_fdd, "dirichlet");

    auto AddBoundarySegment
        = [&mesh] (netgen::PointIndex p0, netgen::PointIndex p1,
                   int edge_descriptor_index) {
            netgen::Segment segment;
            segment[0] = p0;
            segment[1] = p1;
            segment.SetIndex (edge_descriptor_index);
            mesh->AddSegment (segment);
          };

    for (int i = 0; i < tps->GetNTents (); i++)
      {
        const Tent *tent = &tps->GetTent (i);

        if (tent->ttop < tent->tbot + dt_eps)
          {
            cout << "had to skip degenerate tent" << endl;
            continue;
          }

        // Add vertices and 2d elements to the mesh
        Vector<netgen::PointIndex> vertices (tent->nbv.Size () + 2);
        const double pointc = ma->GetPoint<1> (tent->vertex) (0);
        const int ibot = tent->nbv[0] > tent->vertex ? 0 : 2;
        const int itop = tent->nbv[0] > tent->vertex ? 2 : 0;

        vertices[ibot] = AddPointUnique (
            mesh, &pim, netgen::Point3d (pointc, tent->tbot, 0));
        vertices[itop] = AddPointUnique (
            mesh, &pim, netgen::Point3d (pointc, tent->ttop, 0));

        for (size_t k = 0; k < tent->nbv.Size (); k++)
          vertices[2 * k + 1] = AddPointUnique (
              mesh, &pim,
              netgen::Point3d (ma->GetPoint<1> (tent->nbv[k]) (0),
                               tent->nbtime[k], 0));

        if (tent->nbv.Size () == 1)
          {
            netgen::Element2d element (netgen::TRIG);
            for (int j = 0; j < 3; j++)
              element[j] = vertices[j];
            element.SetIndex (domain_index);
            mesh->AddSurfaceElement (element);
          }
        else if (tent->nbv.Size () == 2)
          {
            netgen::Element2d element (netgen::QUAD);
            for (int j = 0; j < 4; j++)
              element[j] = vertices[j];
            element.SetIndex (domain_index);
            mesh->AddSurfaceElement (element);
          }

        for (size_t k = 0; k < tent->nbv.Size (); k++)
          {
            // Add 1d Elements - inflow
            if (tent->tbot < dt_eps && tent->nbtime[k] < dt_eps)
              {
                if (tent->vertex < tent->nbv[k])
                  AddBoundarySegment (vertices[ibot], vertices[2 * k + 1],
                                      ed_inflow);
                else
                  AddBoundarySegment (vertices[2 * k + 1], vertices[ibot],
                                      ed_inflow);
              }
            // Add 1d Elements - outflow
            if (tent->ttop > dt - dt_eps && tent->nbtime[k] > dt - dt_eps)
              {
                if (tent->vertex < tent->nbv[k])
                  AddBoundarySegment (vertices[2 * k + 1], vertices[itop],
                                      ed_outflow);
                else
                  AddBoundarySegment (vertices[itop], vertices[2 * k + 1],
                                      ed_outflow);
              }
            // Add 1d Elements - dirichlet
            if (pointc == bd_points[0] || pointc == bd_points[1])
              AddBoundarySegment (vertices[2], vertices[0], ed_dirichlet);
          }
      }

    mesh->SetBCName (0, "inflow");
    mesh->SetBCName (1, "outflow");
    mesh->SetBCName (2, "dirichlet");

    return make_shared<MeshAccess> (mesh);
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
