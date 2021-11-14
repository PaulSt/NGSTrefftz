#include <python_ngstd.hpp>
#include "python_fem.hpp"
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
