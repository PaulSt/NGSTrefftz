#ifndef FILE_HELPERS_HPP
#define FILE_HELPERS_HPP

#include <cmath>
#include <bla.hpp>
#include <fem.hpp>

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

  template <int D, typename SCAL>
  inline SCAL Dot (const AutoDiff<D, SCAL> &u, const AutoDiff<D, SCAL> &v)
  {
    SCAL sum = 0;
    for (int i = 0; i < D; i++)
      sum += u.DValue (i) * v.DValue (i);
    return sum;
  }
}

#endif
