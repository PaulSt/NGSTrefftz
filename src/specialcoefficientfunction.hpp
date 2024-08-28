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

  /**
   * \class PrintCF
   * \brief A class that represents a coefficient function for printing
   * integration point values to a CSV file.
   *
   * This class inherits from the CoefficientFunction class and provides
   * functionality for writing the the world coordinates of evaluated mapped
   * integration points to a CSV file.
   */
  class PrintCF : public CoefficientFunction
  {
    string filename = "printcf.csv";
    shared_ptr<ofstream> ofs = nullptr;

  public:
    PrintCF (const string &a_filename)
    {
      filename = a_filename;
      ofs = make_shared<ofstream> (filename);
    }

    double Evaluate (const BaseMappedIntegrationPoint &mip) const;

    ~PrintCF () { ofs->close (); }
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSpecialCoefficientFunction (py::module m);
#endif // NGS_PYTHON

#endif
