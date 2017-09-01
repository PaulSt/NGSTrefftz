#ifndef __FILE_TREFFTZCOEFFICIENT_HPP
#define __FILE_TREFFTZCOEFFICIENT_HPP
#include "TrefftzElement.hpp"

namespace ngfem
{
  class TrefftzCoefficientFunction : public CoefficientFunction
  {
    int basisfunction;

  public:
    TrefftzCoefficientFunction () : CoefficientFunction (1) { ; }

    TrefftzCoefficientFunction (int basis) : CoefficientFunction (1)
    {
      basisfunction = basis;
    }

    virtual double
    Evaluate (const BaseMappedIntegrationPoint &mip) const override
    {
      FlatVector<double> point = mip.GetPoint ();

      TrefftzElement<2, 4> treff;

      int ndof = treff.GetNBasis ();

      Vector<> shape (ndof);

      treff.CalcShape (mip, shape);

      return shape (basisfunction); // shape(0);
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzCoefficient (py::module m)
{
  py::class_<TrefftzCoefficientFunction,
             shared_ptr<TrefftzCoefficientFunction>, CoefficientFunction> (
      m, "TrefftzCoefficient", "")
      .def (py::init<> ())
      .def (py::init<int> ());
}
#endif // NGS_PYTHON

#endif //  __FILE_MYCOEFFICIENT_HPP
