#ifndef __FILE_TREFFTZCOEFFICIENT_HPP
#define __FILE_TREFFTZCOEFFICIENT_HPP
#include "TrefftzElement.hpp"

namespace ngfem
{
  class TrefftzCoefficientFunction : public CoefficientFunction
  {
		int basisfunction;
		TrefftzElement<2,3> treff;

  public:
    TrefftzCoefficientFunction()
      : CoefficientFunction(1) { ; }

		TrefftzCoefficientFunction(int basis)
			: CoefficientFunction(1) { basisfunction = basis; }

    virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
    {
			FlatVector<double> point = mip.GetPoint();

			int ndof = treff.GetNBasis();

			Vector<double> shape(ndof);
			//Matrix<double> shape(ndof,2);
			treff.CalcShape(mip,shape);
      return shape(basisfunction);
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzCoefficient(py::module m)
{
  py::class_<TrefftzCoefficientFunction, shared_ptr<TrefftzCoefficientFunction>, CoefficientFunction>
    (m, "TrefftzCoefficient", "")
    .def(py::init<>())
		.def(py::init<int>())
    ;
}
#endif // NGS_PYTHON

#endif //  __FILE_MYCOEFFICIENT_HPP
