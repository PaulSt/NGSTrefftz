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
  int EvolveTents (shared_ptr<MeshAccess> ma)
  {
    py::module et = py::module::import ("EvolveTent");
    int basedim = ma->GetDimension ();
    TentPitchedSlab<1> tps
        = TentPitchedSlab<1> (ma); // collection of tents in timeslab
    if (basedim >= 2)
      cout << "oh no" << endl;
    tps.PitchTents (1.0, 1.0); // adt = time slab height, wavespeed

    auto mesh = make_shared<netgen::Mesh> ();
    mesh->SetDimension (ma->GetDimension () + 1);
    netgen::SetGlobalMesh (mesh); // for visualization
    mesh->SetGeometry (make_shared<netgen::NetgenGeometry> ());
    netgen::PointIndex pind;

    for (Tent *tent : tps.tents)
      {
        pind = mesh->AddPoint (netgen::Point3d (1, 1, 1));
        // cout << *tent << endl;
        for (auto v : tent->nbv)
          ma->GetPoint<1> (v);
      }
    cout << endl;

    // Flags flag();
    // flag.SetFlag("order",4);
    TrefftzFESpace fes (ma, Flags ());
    // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
    // py::object ffes = py::cast(fes);
    // auto pyspace = py::class_<TrefftzFESpace,
    // shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str());
    py::object pyfes = et.attr ("GetFESTrefftz") (ma);
    // FESpace *ffes = pyfes.cast<FESpace *>();
    // et.attr("EvolveTent")(pyfes,?,?);
    return 4;
  }
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m)
{
  m.def ("EvolveTents",
         [] (shared_ptr<MeshAccess> ma) -> int {
           return EvolveTents (ma);
         } //, py::call_guard<py::gil_scoped_release>()
  );
}
#endif // NGS_PYTHON
