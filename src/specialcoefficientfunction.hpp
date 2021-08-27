#ifndef SPECIALCOEFFICIENTFUNCTION_HPP
#define SPECIALCOEFFICIENTFUNCTION_HPP

#include <fem.hpp>
#include <comp.hpp>
#include <multigrid.hpp>
#include <h1lofe.hpp>
#include <regex>

using namespace ngcomp;
namespace ngfem
{

    class ClipCoefficientFunction : public CoefficientFunction
    {
        private:
            shared_ptr<CoefficientFunction> coef;
            int clipdim;
            double clipvalue;
        public:
            ClipCoefficientFunction(shared_ptr<CoefficientFunction> acoef,int adimension,int aclipdim,double aclipvalue, bool ais_complex = false)
                : CoefficientFunction(adimension,ais_complex), coef(acoef), clipdim(aclipdim), clipvalue(aclipvalue)
            { ; }
            using CoefficientFunction::Evaluate;
            virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
            virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
            virtual void EvaluateStdRule (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    };


    class IntegrationPointFunction : public CoefficientFunction
    {
        private:
            vector<vector<double>> values;
        public:
            IntegrationPointFunction(shared_ptr<MeshAccess> mesh, IntegrationRule& intrule, Vector<> ipdata);
            IntegrationPointFunction(shared_ptr<MeshAccess> mesh, IntegrationRule& intrule, Matrix<> ipdata);
            virtual double Evaluate(const BaseMappedIntegrationPoint & ip) const;
            void PrintTable();
    };


    class WeightedRadiusFunction : public CoefficientFunction
    {
        private:
            Vector<> values;
        public:
            WeightedRadiusFunction(shared_ptr<MeshAccess> mesh, shared_ptr<CoefficientFunction> wavespeedcf);
            virtual double Evaluate(const BaseMappedIntegrationPoint & ip) const;
    };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSpecialCoefficientFunction(py::module m);
#endif // NGS_PYTHON

#endif
