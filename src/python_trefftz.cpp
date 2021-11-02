#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
using namespace ngsolve;
#include "python_fem.hpp"
#include "trefftzfespace.hpp"
#include "monomialfespace.hpp"
#include "diffopmapped.hpp"
#include "specialcoefficientfunction.hpp"
#include <tents.hpp>
#include "meshtentslab.hpp"
#include "twavetents.hpp"
#include <boost/math/special_functions/airy.hpp>

namespace ngcomp
{
  // Airy function used to test quasi-Trefftz methods
  struct GenericAiry
  {
    double operator() (double x) const { return boost::math::airy_ai (x); }
    SIMD<double> operator() (SIMD<double> x) const
    {
      return SIMD<double> (
          [&] (int i) -> double { return boost::math::airy_ai (x[i]); });
    }
    template <typename T> T operator() (T x) const
    {
      throw Exception (string ("airy not available for type ")
                       + typeid (T).name ());
    }
    template <typename T> AutoDiff<1, T> operator() (AutoDiff<1, T> x) const
    {
      throw Exception ("airy(..) is not complex differentiable");
    }
    template <typename T>
    AutoDiffDiff<1, T> operator() (AutoDiffDiff<1, T> x) const
    {
      throw Exception ("airy(..) is not complex differentiable");
    }
    static string Name () { return "airy"; }
    void DoArchive (Archive &ar) {}
  };

  struct GenericAiryP
  {
    double operator() (double x) const
    {
      return boost::math::airy_ai_prime (x);
    }
    SIMD<double> operator() (SIMD<double> x) const
    {
      return SIMD<double> (
          [&] (int i) -> double { return boost::math::airy_ai_prime (x[i]); });
    }
    template <typename T> T operator() (T x) const
    {
      throw Exception (string ("airy prime not available for type ")
                       + typeid (T).name ());
    }
    template <typename T> AutoDiff<1, T> operator() (AutoDiff<1, T> x) const
    {
      throw Exception ("airyp(..) is not complex differentiable");
    }
    template <typename T>
    AutoDiffDiff<1, T> operator() (AutoDiffDiff<1, T> x) const
    {
      throw Exception ("airyp(..) is not complex differentiable");
    }
    static string Name () { return "airyp"; }
    void DoArchive (Archive &ar) {}
  };
}

PYBIND11_MODULE (tngs, m)
{
  py::module::import ("ngsolve");
  // py::module::import("ngstents");
  m.attr ("__name__") = "tngs";
  m.attr ("__package__") = "tngs";

  ExportTrefftzFESpace (m);
  ExportMonomialFESpace (m);
  ExportSpecialCoefficientFunction (m);
  ExportTWaveTents (m);
  ExportMeshTentSlab (m);
  ExportStdMathFunction<GenericAiry> (m, "airy", "airy function");
  ExportStdMathFunction<GenericAiryP> (m, "airyp", "airyp function");
}
