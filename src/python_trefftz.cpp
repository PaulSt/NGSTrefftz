#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
#include <comp.hpp>
// using namespace ngsolve;
#include "trefftzfespace.hpp"
#include "specialcoefficientfunction.hpp"
#include "specialintegrator.hpp"
#ifdef NGSTREFFTZ_USE_NGSTENTS
#include "twavetents.hpp"
#include "mesh1dtents.hpp"
#endif
#include "embtrefftz.hpp"
#include "monomialfespace.hpp"
#include "pufespace.hpp"
#include "condensedg.hpp"
#include "boxintegral.hpp"
#include "tp0fespace.hpp"
// #include "airy.cpp"

PYBIND11_MODULE (ngstrefftz, m)
{
  py::module::import ("ngsolve");
  m.attr ("__name__") = "ngstrefftz";
  m.attr ("__package__") = "ngstrefftz";

  ExportTrefftzFESpace (m);
  ExportSpecialCoefficientFunction (m);
  ExportSpecialIntegrator (m);
#ifdef NGSTREFFTZ_USE_NGSTENTS
  ExportTWaveTents (m);
  ExportMesh1dTents (m);
#endif
  ExportEmbTrefftz (m);
  ExportMonomialFESpace (m);
  ExportPUFESpace (m);
  ExportCondenseDG (m);
  ExportBoxIntegral (m);
  ExportTP0FESpace (m);
  // ExportStdMathFunction<GenericAiry>(m, "airy", "airy function");
  // ExportStdMathFunction<GenericAiryP>(m, "airyp", "airyp function");
}
