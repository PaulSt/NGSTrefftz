#include <solve.hpp>
using namespace ngsolve;
#include <python_ngstd.hpp>
//#include "trefftz/trefftzelement.hpp"
#include "trefftz/trefftzfespace.hpp"
#include "trefftz/trefftzcoefficient.hpp"
#include "trefftz/diffopmapped.hpp"
#include "trefftz/clipcoefficientfunction.hpp"
//#include "trefftz/mappedelement.hpp"
#include "trefftz/tents/tents.hpp"
#include "trefftz/evolvetent.hpp"
#include "trefftz/meshtentslab.hpp"

PYBIND11_PLUGIN(trefftzngs) {
    // import ngsolve such that python base classes are defined
    py::module::import("ngsolve");

    py::module m("trefftzngs", "trefftzngs doc-string");

    //ExportMappedElement(m);
    ExportTrefftzFESpace(m);
    ExportTrefftzCoefficient(m);
    ExportClipCoefficientFunction(m);
    ExportEvolveTent(m);
    ExportMeshTentSlab(m);

    return m.ptr();
}

//static RegisterNumProc<NumProcPyDemo> npinit1("demopy");
