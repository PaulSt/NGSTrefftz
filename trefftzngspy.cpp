#include <solve.hpp>
using namespace ngsolve;
#include <python_ngstd.hpp>
#include "trefftz/TrefftzElement.hpp"
#include "trefftz/TrefftzFESpace.hpp"
#include "trefftz/MultiArray.hpp"
#include "trefftz/TrefftzCoefficient.hpp"

PYBIND11_PLUGIN (trefftzngs)
{
  // import ngsolve such that python base classes are defined
  py::module::import ("ngsolve");

  py::module m ("trefftzngs", "trefftzngs doc-string");

  /*
    Finally we export all the other classes and functions we created in this
    tutorial
   */
  ExportTrefftzElement (m);
  ExportMultiArray (m);
  ExportTrefftzFESpace (m);
  ExportTrefftzCoefficient (m);

  return m.ptr ();
}

// static RegisterNumProc<NumProcPyDemo> npinit1("demopy");
