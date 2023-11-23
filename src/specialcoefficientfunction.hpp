#ifndef SPECIALCOEFFICIENTFUNCTION_HPP
#define SPECIALCOEFFICIENTFUNCTION_HPP

#include <fem.hpp>
#include <comp.hpp>
#include <multigrid.hpp>
#include <h1lofe.hpp>
#include <regex>

namespace ngfem
{
  using namespace ngcomp;

  class ClipCoefficientFunction : public CoefficientFunction
  {
  private:
    shared_ptr<CoefficientFunction> coef;
    int clipdim;
    double clipvalue;

  public:
    ClipCoefficientFunction (shared_ptr<CoefficientFunction> acoef,
                             int adimension, int aclipdim, double aclipvalue,
                             bool ais_complex = false)
        : CoefficientFunction (adimension, ais_complex), coef (acoef),
          clipdim (aclipdim), clipvalue (aclipvalue)
    {
      ;
    }
    using CoefficientFunction::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint &ip) const override;
    void Evaluate (const BaseMappedIntegrationPoint &ip,
                   FlatVector<> result) const override;
  };

  class IntegrationPointFunction : public CoefficientFunction
  {
  private:
    vector<vector<double>> values;
    shared_ptr<MeshAccess> ma;
    IntegrationRule intrule;

  public:
    IntegrationPointFunction (shared_ptr<MeshAccess> mesh,
                              IntegrationRule &intrule, Vector<> ipdata);
    IntegrationPointFunction (shared_ptr<MeshAccess> mesh,
                              IntegrationRule &intrule, Matrix<> ipdata);
    virtual double Evaluate (const BaseMappedIntegrationPoint &ip) const;
    void PrintTable ();
    vector<vector<double>> Export ();
  };

  class WeightedRadiusFunction : public CoefficientFunction
  {
  private:
    Vector<> values;

  public:
    WeightedRadiusFunction (shared_ptr<MeshAccess> mesh,
                            shared_ptr<CoefficientFunction> wavespeedcf);
    virtual double Evaluate (const BaseMappedIntegrationPoint &ip) const;
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSpecialCoefficientFunction (py::module m);
#endif // NGS_PYTHON

#endif
