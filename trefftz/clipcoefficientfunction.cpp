#include <fem.hpp>
#include "clipcoefficientfunction.hpp"
#include <comp.hpp>
using namespace ngcomp;

namespace ngfem
{

double ClipCoefficientFunction :: Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
		ip.GetPoint()[clipdim] = clipvalue;
    return coef->Evaluate(ip);
  }

  void ClipCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
		cout << "not yet implemented" << endl;
  }

  void ClipCoefficientFunction :: EvaluateStdRule (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    for(int i=0;i<ir.Size();i++)
      values(i,0) = Evaluate(ir[i]);
  }
}


typedef CoefficientFunction CF;
typedef shared_ptr<CoefficientFunction> spCF;

#ifdef NGS_PYTHON
void ExportClipCoefficientFunction(py::module m)
{
	m.def("ClipCoefficientFunction", [](spCF cf_x, int aclipdim,double aclipvalue) -> shared_ptr<CoefficientFunction>
					{
						auto pcf = make_shared<ClipCoefficientFunction>(cf_x,cf_x->Dimension(),aclipdim,aclipvalue,false);
						pcf->SetDimension(pcf->Dimension());
						return pcf;
					},
			 py::call_guard<py::gil_scoped_release>());
}
#endif // NGS_PYTHON
