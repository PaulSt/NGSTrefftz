#ifndef FILE_HELPERS_HPP
#define FILE_HELPERS_HPP

#include <cmath>
#include <bla.hpp>
#include <fem.hpp>
#include <boost/math/special_functions/airy.hpp>

namespace ngfem
{
  using namespace ngbla;

  typedef Vec<3, Array<double>>
      CSR; // CSR sparse matrix in (row,col,val) format

  void MatToCSR (Matrix<> mat, CSR &sparsemat);

  class Monomial : public RecursivePolynomial<Monomial>
  {
  public:
    Monomial () { ; }

    template <class S, class T> inline Monomial (int n, S x, T &&values)
    {
      Eval (n, x, values);
    }

    template <class S> static INLINE double P0 (S x) { return 1.0; }
    template <class S> static INLINE S P1 (S x) { return x; }
    template <class S, class Sy> static INLINE S P1 (S x, Sy y)
    {
      return P1 (x);
    }

    static INLINE double A (int i) { return 1.0; }
    static INLINE double B (int i) { return 0; }
    static INLINE double C (int i) { return 0; }

    static INLINE double CalcA (int i) { return 1.0; }
    static INLINE double CalcB (int i) { return 0; }
    static INLINE double CalcC (int i) { return 0; }

    enum
    {
      ZERO_B = 1
    };
  };

  constexpr long ipow (int base, int expo, int result = 1)
  {
    return expo < 1 ? result
                    : ipow (base * base, expo / 2,
                            (expo % 2) ? result * base : result);
  }

  inline constexpr long BinCoeff (int n, int k)
  {
    return round (tgamma (n + 1) / (tgamma (k + 1) * tgamma (n - k + 1)));
  }

  // k-th coeff of Legendre polynomial of degree n in monomial basis
  constexpr double LegCoeffMonBasis (int n, int k)
  {
    if (n == 0)
      return 1;
    if (k > n)
      return 0;
    if ((n + k) % 2)
      return 0;
    double coeff = pow (2, -n) * pow (-1, floor ((n - k) / 2))
                   * BinCoeff (n, floor ((n - k) / 2)) * BinCoeff (n + k, n);
    // double coeff = pow(2,-n) * pow(-1,k) * BinCoeff(n,k) *
    // BinCoeff(2*n-2*k,n);
    return coeff;
  }

  // k-th coeff of Chebyshev polynomial of degree n in monomial basis
  constexpr double ChebCoeffMonBasis (int n, int k)
  {
    if (n == 0)
      return 1;
    if (k > n)
      return 0;
    if ((n + k) % 2)
      return 0;
    double coeff = pow (2, k - 1) * n * pow (-1, floor ((n - k) / 2))
                   * tgamma ((n + k) / 2)
                   / (tgamma (floor ((n - k) / 2) + 1) * tgamma (k + 1));
    return coeff;
  }

  template <typename T> int sgn_nozero (T val)
  {
    return (T (0) <= val) - (val < T (0));
  }

  constexpr int factorial (int n) { return n > 1 ? n * factorial (n - 1) : 1; }

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

#endif
