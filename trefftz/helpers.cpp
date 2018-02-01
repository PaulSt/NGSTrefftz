#ifndef FILE_HELPERS
#define FILE_HELPERS
#include <cmath>

namespace ngfem
{
	constexpr long ipow(int base, int expo, int result = 1)
	{
		return expo < 1 ? result : ipow(base*base, expo/2, (expo % 2) ? result*base : result);
	}

	inline constexpr long BinCoeff(int n,int k)
	{
		return round( tgamma(n+1) / (tgamma(k+1)*tgamma(n-k+1)) );
	}

	// k-th coeff of Legendre polynomial of degree n in monomial basis
	constexpr double LegCoeffMonBasis(int n, int k)
	{
		if(k>n) return 0;
		if((n+k)%2) return 0;
		double coeff = pow(2,-n) * pow(-1,floor((n-k)/2)) * BinCoeff(n,floor((n-k)/2)) * BinCoeff(n+k,n);
		// double coeff = pow(2,-n) * pow(-1,k) * BinCoeff(n,k) * BinCoeff(2*n-2*k,n);
		return coeff;
	}
}

#endif
