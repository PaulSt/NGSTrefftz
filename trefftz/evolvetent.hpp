#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <Python.h>
#include <python_ngstd.hpp>
#include <pybind11/iostream.h>

namespace ngcomp
{
  int TestPython (shared_ptr<MeshAccess> ma)
  {
    py::module et = py::module::import ("EvolveTent");
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
void ExportEvolveTent (py::module m)
{
  m.def ("TestPythonF",
         [] (shared_ptr<MeshAccess> ma) -> int {
           return TestPython (ma);
         } //, py::call_guard<py::gil_scoped_release>()
  );
};
#endif // NGS_PYTHON

#endif
