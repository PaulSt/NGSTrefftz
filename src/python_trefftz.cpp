#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
#include <comp.hpp>
// using namespace ngsolve;
#include <tents.hpp>
#include <python_tents.cpp>
#include "trefftzfespace.hpp"
#include "specialcoefficientfunction.hpp"
#include "specialintegrator.hpp"
#include "twavetents.hpp"
#include "embtrefftz.hpp"
#include "mesh1dtents.hpp"
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

  ExportTents (m);
  ExportTrefftzFESpace (m);
  ExportSpecialCoefficientFunction (m);
  ExportSpecialIntegrator (m);
  ExportTWaveTents (m);
  ExportEmbTrefftz (m);
  ExportMesh1dTents (m);
  ExportMonomialFESpace (m);
  ExportPUFESpace (m);
  ExportCondenseDG (m);
  ExportBoxIntegral (m);
  ExportTP0FESpace (m);
  // ExportStdMathFunction<GenericAiry>(m, "airy", "airy function");
  // ExportStdMathFunction<GenericAiryP>(m, "airyp", "airyp function");
}
