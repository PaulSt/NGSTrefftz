#include "evolvetent.hpp"
#include "trefftzelement.hpp"
#include "tents/tents.hpp"
#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>

namespace ngcomp
{
  template <int D>
  void EvolveTents (shared_ptr<MeshAccess> ma, double wavespeed, double dt,
                    Vector<double> wavefront)
  {
    int order = 3;
    int ir_order = ceil ((order + 1) / 1);
    int ne = ma->GetNE ();

    // wavefront.Resize(ir_order * ma->GetNEdges() * (D+2));
    T_TrefftzElement<D + 1> tel (order, wavespeed);
    IntegrationRule ir (ET_SEGM, ir_order);
    ScalarFE<ET_SEGM, D> faceint;

    // py::module et = py::module::import("DGeq");

    TentPitchedSlab<D> tps
        = TentPitchedSlab<1> (ma);  // collection of tents in timeslab
    tps.PitchTents (dt, wavespeed); // adt = time slab height, wavespeed
    LocalHeap lh (order * D * 1000);

    RunParallelDependency (tps.tent_dependency, [&] (int i) {
      HeapReset hr (lh);
      // LocalHeap slh = lh.Split();  // split to threads
      Tent &tent = *tps.tents[i];
      // cout << tent << endl;
      cout << endl << "tent: " << i << " " << tent.vertex << endl;

      int nbasis = BinCoeff (D - 1 + order, order)
                   + BinCoeff (D - 1 + order - 1, order - 1);
      Matrix<double> elmatrix (nbasis, nbasis);
      Vector<double> elvector (nbasis);
      // for (int k = 0; k < tent.nbv.Size(); k++)
      for (auto el : tent.els)
        {
          auto verts = ma->GetEdgePNums (el);
          MappedIntegrationRule<1, D> mir (ir, ma->GetTrafo (el, lh),
                                           lh); // <dim  el, dim space>

          // determine linear basis function to use for tent face
          Vec<3> bs = 0;
          for (int ivert = 0; ivert < verts.Size (); ivert++)
            {
              if (verts[ivert] == tent.vertex)
                bs[ivert] = 1;
            }
          auto p = ma->GetPoint<1> (tent.vertex);
          cout << "mir0: " << mir[0]
               << " val: " << faceint.Evaluate (ir[0], bs) << endl;
          cout << "mir1: " << mir[1]
               << " val: " << faceint.Evaluate (ir[1], bs) << endl;
          // cout << "mir1: " << mir[ir.Size()-1]<<mir[ir.Size()-1].GetWeight()
          // << endl;

          for (int imip = 0; imip < mir.Size (); imip++)
            {
              auto p = mir[imip].GetPoint ();
              Vec<D + 1> tmip;
              tmip.Range (0, D) = p;
            }
        }
    });
    // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
    // py::object ffes = py::cast(fes);
    // auto pyspace = py::class_<TrefftzFESpace,
    // shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str()); py::object
    // pyfes = et.attr("GetFESTrefftz")(ma); FESpace *ffes = pyfes.cast<FESpace
    // *>(); et.attr("EvolveTent")(pyfes,?,?);
  }

  shared_ptr<MeshAccess>
  NgsTPmesh (shared_ptr<MeshAccess> ma, double wavespeed, double dt)
  {
    Point2IndexMap *pim = new Point2IndexMap (); // Your map type may vary,
                                                 // just change the typedef
    int index = 1;
    int basedim = ma->GetDimension ();
    double dt_eps = 0.00001;
    // get boundaries of the init mesh
    Array<double> bd_points (0);
    for (auto el : ma->Elements (BND))
      bd_points.Append (ma->GetPoint<1> (el.Vertices ()[0]) (0));

    TentPitchedSlab<1> tps
        = TentPitchedSlab<1> (ma);  // collection of tents in timeslab
    tps.PitchTents (dt, wavespeed); // adt = time slab height, wavespeed
    cout << "numb of tents: " << tps.tents.Size () << endl;

    auto mesh = make_shared<netgen::Mesh> ();
    mesh->SetDimension (ma->GetDimension () + 1);
    netgen::SetGlobalMesh (mesh); // for visualization
    mesh->SetGeometry (make_shared<netgen::NetgenGeometry> ());
    mesh->SetMaterial (1, "mat");

    netgen::FaceDescriptor fd (1, 1, 0, 1);
    int ind_fd = mesh->AddFaceDescriptor (fd);

    netgen::FaceDescriptor fdi (1, 1, 0, 1);
    fdi.SetBCName (new string ("inflow"));
    int ind_fdi = mesh->AddFaceDescriptor (fdi);

    netgen::FaceDescriptor fdo (2, 1, 0, 2);
    fdo.SetBCName (new string ("outflow"));
    int ind_fdo = mesh->AddFaceDescriptor (fdo);

    netgen::FaceDescriptor fdd (3, 1, 0, 3);
    fdd.SetBCName (new string ("dirichlet"));
    int ind_fdd = mesh->AddFaceDescriptor (fdd);

    for (Tent *tent : tps.tents)
      {
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
        for (int k = 0; k < tent->nbv.Size (); k++)
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

        for (int k = 0; k < tent->nbv.Size (); k++)
          {
            // Add 1d Elements - inflow
            if (tent->tbot < dt_eps && tent->nbtime[k] < dt_eps)
              {
                netgen::Segment *newel = new netgen::Segment ();
                (*newel)[tent->vertex < tent->nbv[k] ? 0 : 1] = vertices[ibot];
                (*newel)[tent->vertex < tent->nbv[k] ? 1 : 0]
                    = vertices[2 * k + 1];
                newel->si = ind_fdi;
                newel->edgenr = 1;
                newel->epgeominfo[0].edgenr = 1;
                newel->epgeominfo[1].edgenr = 1;
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
    if (pi == pim->size ())
      ngma->AddPoint (p);
    return pi;
  }

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m)
{
  m.def ("EvolveTents",
         [] (shared_ptr<MeshAccess> ma, double wavespeed,
             double dt) //-> shared_ptr<MeshAccess>
         {
           int D = 1;
           int order = 3;
           int ir_order = ceil ((order + 1) / 1);
           int ne = ma->GetNE ();
           Vector<double> wavefront (ir_order * ne * (D + 2));
           return EvolveTents<1> (ma, wavespeed, dt, wavefront);
         } //, py::call_guard<py::gil_scoped_release>()
  );
  m.def ("NgsTPmesh",
         [] (shared_ptr<MeshAccess> ma, double wavespeed,
             double dt) -> shared_ptr<MeshAccess> {
           return NgsTPmesh (ma, wavespeed, dt);
         } //, py::call_guard<py::gil_scoped_release>()
  );
}
#endif // NGS_PYTHON
