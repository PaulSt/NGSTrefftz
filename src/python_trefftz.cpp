#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
using namespace ngsolve;
#include "python_fem.hpp"
#include "trefftzfespace.hpp"
#include "monomialfespace.hpp"
#include "diffopmapped.hpp"
#include "specialcoefficientfunction.hpp"
#include <tents.hpp>
#include "meshtentslab.hpp"
#include "evolvetent.hpp"

PYBIND11_MODULE(trefftzngs,m) {
    py::module::import("ngsolve");
    //py::module::import("ngstents");
    //m.attr("__name__") = "trefftzngs"
    m.attr("__package__") = "trefftzngs";

    ExportTrefftzFESpace(m);
    ExportMonomialFESpace(m);
    ExportSpecialCoefficientFunction(m);
    ExportEvolveTent(m);
    ExportMeshTentSlab(m);
    ExportStdMathFunction<GenericAiry>(m, "airy", "airy function");
    ExportStdMathFunction<GenericAiryP>(m, "airyp", "airyp function");
}

