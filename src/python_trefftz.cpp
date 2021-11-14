#include <python_ngstd.hpp>
#include <solve.hpp>
//#include <fem.hpp>
// using namespace ngsolve;
#include "trefftzfespace.hpp"
#include "monomialfespace.hpp"
#include "diffopmapped.hpp"
#include "specialcoefficientfunction.hpp"
//#include <tents.hpp>
#include "meshtentslab.hpp"
#include "twavetents.hpp"
//#include "svdtrefftz.hpp"

#include "airy.cpp"

PYBIND11_MODULE (_tngs, m)
{
  py::module::import ("ngsolve");
  m.attr ("__name__") = "tngs";
  m.attr ("__package__") = "tngs";

  ExportTrefftzFESpace (m);
  ExportMonomialFESpace (m);
  ExportSpecialCoefficientFunction (m);
  ExportTWaveTents (m);
  ExportMeshTentSlab (m);

  ExportStdMathFunction<GenericAiry> (m, "airy", "airy function");
  ExportStdMathFunction<GenericAiryP> (m, "airyp", "airyp function");
}
