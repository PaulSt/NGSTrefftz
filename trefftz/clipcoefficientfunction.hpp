#ifndef CLIPCOEFFICIENTFUNCTION_HPP
#define CLIPCOEFFICIENTFUNCTION_HPP

namespace ngfem
{

  class ClipCoefficientFunction : public CoefficientFunction
  {
  private:
    shared_ptr<CoefficientFunction> coef;
    double clipvalue;
    int clipdim;

  public:
    ///
    ClipCoefficientFunction (shared_ptr<CoefficientFunction> acoef,
                             int adimension, int aclipdim, double aclipvalue,
                             bool ais_complex = false)
        : CoefficientFunction (adimension, ais_complex), coef (acoef),
          clipdim (aclipdim), clipvalue (aclipvalue)
    {
      ;
    }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint &ip) const;
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule &ir,
                           FlatMatrix<double> values) const;
    virtual void EvaluateStdRule (const BaseMappedIntegrationRule &ir,
                                  FlatMatrix<double> values) const;
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportClipCoefficientFunction (py::module m);
#endif // NGS_PYTHON

#endif
