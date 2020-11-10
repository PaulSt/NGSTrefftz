#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
using namespace ngsolve;
#include "python_fem.hpp"
//#include "trefftz/trefftzwavefe.hpp"
#include "trefftz/trefftzfespace.hpp"
#include "trefftz/diffopmapped.hpp"

PYBIND11_PLUGIN(trefftzngs) {
    // import ngsolve such that python base classes are defined
    py::module::import("ngsolve");

    py::module m("trefftzngs", "trefftzngs doc-string");

    ExportTrefftzFESpace(m);
    ExportStdMathFunction<GenericAiry>(m, "airy", "airy function");
    ExportStdMathFunction<GenericAiryP>(m, "airyp", "airyp function");

    return m.ptr();
}

//static RegisterNumProc<NumProcPyDemo> npinit1("demopy");
