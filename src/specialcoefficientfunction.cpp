#include "specialcoefficientfunction.hpp"
#include <fem.hpp>
#include <comp.hpp>
#include <h1lofe.hpp>
#include <regex>
#include <multigrid.hpp>
using namespace ngcomp;

namespace ngfem
{

  double ClipCoefficientFunction ::Evaluate (
      const BaseMappedIntegrationPoint &ip) const
  {
    ip.GetPoint ()[clipdim] = clipvalue;
    return coef->Evaluate (ip);
  }

  void ClipCoefficientFunction ::Evaluate (const BaseMappedIntegrationRule &ir,
                                           FlatMatrix<double> values) const
  {
    cout << "not yet implemented" << endl;
  }

  void ClipCoefficientFunction ::EvaluateStdRule (
      const BaseMappedIntegrationRule &ir, FlatMatrix<double> values) const
  {
    for (int i = 0; i < ir.Size (); i++)
      values (i, 0) = Evaluate (ir[i]);
  }
}

typedef CoefficientFunction CF;
typedef shared_ptr<CoefficientFunction> spCF;

#ifdef NGS_PYTHON
void ExportSpecialCoefficientFunction (py::module m)
{
  m.def (
      "ClipCoefficientFunction",
      [] (spCF cf_x, int aclipdim,
          double aclipvalue) -> shared_ptr<CoefficientFunction> {
        auto pcf = make_shared<ClipCoefficientFunction> (
            cf_x, cf_x->Dimension (), aclipdim, aclipvalue, false);
        pcf->SetDimension (pcf->Dimension ());
        return pcf;
      },
      py::call_guard<py::gil_scoped_release> ());

  py::class_<IntegrationPointFunction, shared_ptr<IntegrationPointFunction>,
             CoefficientFunction> (m, "IntegrationPointFunction")
      .def (
          "__init__",
          [] (IntegrationPointFunction *instance, shared_ptr<MeshAccess> mesh,
              IntegrationRule &intrule, Vector<> data) {
            new (instance) IntegrationPointFunction (mesh, intrule, data);
          },
          py::arg ("mesh"), py::arg ("intrule"), py::arg ("Vector"))
      .def (
          "__init__",
          [] (IntegrationPointFunction *instance, shared_ptr<MeshAccess> mesh,
              IntegrationRule &intrule, Matrix<> data) {
            new (instance) IntegrationPointFunction (mesh, intrule, data);
          },
          py::arg ("mesh"), py::arg ("intrule"), py::arg ("Matrix"))
      .def ("PrintTable", &IntegrationPointFunction::PrintTable);

  py::class_<WeightedRadiusFunction, shared_ptr<WeightedRadiusFunction>,
             CoefficientFunction> (m, "WeightedRadiusFunction")
      .def (
          "__init__",
          [] (WeightedRadiusFunction *instance, shared_ptr<MeshAccess> mesh,
              shared_ptr<CoefficientFunction> wavespeedcf) {
            new (instance) WeightedRadiusFunction (mesh, wavespeedcf);
          },
          py::arg ("mesh"), py::arg ("CoefficientFunction"));

  // py::class_<TrefftzCoefficientFunction,
  // shared_ptr<TrefftzCoefficientFunction>, CoefficientFunction> (m,
  //"TrefftzCoefficient", "") .def(py::init<>()) .def(py::init<int>())
  //;
}
#endif // NGS_PYTHON
