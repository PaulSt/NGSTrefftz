#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
#include <comp.hpp>
// using namespace ngsolve;
// #include <tents.hpp>
#include "trefftzfespace.hpp"
#include "specialcoefficientfunction.hpp"
#include "specialintegrator.hpp"
#include "twavetents.hpp"
#include "embtrefftz.hpp"
#include "mesh1dtents.hpp"
#include "monomialfespace.hpp"
#include "condensedg.hpp"
// #include "airy.cpp"

PYBIND11_MODULE (_trefftz, m)
{
  py::module::import ("ngsolve");
  m.attr ("__name__") = "ngstrefftz";
  m.attr ("__package__") = "ngstrefftz";

  ExportTrefftzFESpace (m);
  ExportSpecialCoefficientFunction (m);
  ExportSpecialIntegrator (m);
  ExportTWaveTents (m);
  ExportEmbTrefftz (m);
  ExportMesh1dTents (m);
  ExportMonomialFESpace (m);
  ExportCondenseDG (m);
  // ExportStdMathFunction<GenericAiry>(m, "airy", "airy function");
  // ExportStdMathFunction<GenericAiryP>(m, "airyp", "airyp function");
}
