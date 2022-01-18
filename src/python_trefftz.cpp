#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
#include <comp.hpp>
//using namespace ngsolve;
#include "trefftzfespace.hpp"
#include "monomialfespace.hpp"
#include "diffopmapped.hpp"
#include "specialcoefficientfunction.hpp"
//#include <tents.hpp>
#include "meshtentslab.hpp"
#include "twavetents.hpp"
#include "embtrefftz.hpp"

//#include "airy.cpp"


PYBIND11_MODULE(_trefftz,m) {
    py::module::import("ngsolve");
    m.attr("__name__") = "ngstrefftz";
    m.attr("__package__") = "ngstrefftz";

    ExportTrefftzFESpace(m);
    ExportMonomialFESpace(m);
    ExportSpecialCoefficientFunction(m);
    ExportTWaveTents(m);
    ExportMeshTentSlab(m);
    ExportEmbTrefftz(m);

    //ExportStdMathFunction<GenericAiry>(m, "airy", "airy function");
    //ExportStdMathFunction<GenericAiryP>(m, "airyp", "airyp function");
}

