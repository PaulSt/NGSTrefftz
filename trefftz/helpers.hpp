#ifndef FILE_HELPERS_HPP
#define FILE_HELPERS_HPP

#include <cmath>
#include <bla.hpp>

namespace ngfem
{
	using namespace ngbla;

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
		if(n==0) return 1;
		if(k>n) return 0;
		if((n+k)%2) return 0;
		double coeff = pow(2,-n) * pow(-1,floor((n-k)/2)) * BinCoeff(n,floor((n-k)/2)) * BinCoeff(n+k,n);
		// double coeff = pow(2,-n) * pow(-1,k) * BinCoeff(n,k) * BinCoeff(2*n-2*k,n);
		return coeff;
	}

	// k-th coeff of Chebyshev polynomial of degree n in monomial basis
	constexpr double ChebCoeffMonBasis(int n, int k)
	{
		if(n==0) return 1;
		if(k>n) return 0;
		if((n+k)%2) return 0;
		double coeff = pow(2,k-1)*n*pow(-1,floor((n-k)/2)) * tgamma((n+k)/2)/(tgamma(floor((n-k)/2)+1)*tgamma(k+1));
		return coeff;
	}

	template<int D>
	class HornerScheme
	{
		private:
			int ord;
			int npoly;
			Matrix<int> pascal;
		public:
			HornerScheme(int aord): ord(aord), npoly(BinCoeff(D + ord, ord)), pascal(pascal_sym())
			{;}

			inline double MultiHorner(Vector<double> coeff, Vector<double> &point) const
			{
				for(int j=ord;j>0;j--)
					for(int i=0;i<D;i++)
					 	for(int k=pascal(i+1,j)-1; k>=0; k--){
							coeff( pascal(D+1,j-1)+k ) += point( i ) * coeff( pascal(D+1,j)+pascal(i,j+1)+k );
							coeff( pascal(D+1,j)+pascal(i,j+1)+k ) = 0;
						}
				return coeff(0);
			}

			inline Vector<double> MultiHornerMat(Matrix<double> coeff, Vector<double> &point) const
			{
				for(int j=ord;j>0;j--)
					for(int i=0;i<D;i++)
					 	for(int k=pascal(i+1,j)-1; k>=0; k--){
							coeff.Col( pascal(D+1,j-1)+k ) += point( i ) * coeff.Col( pascal(D+1,j)+pascal(i,j+1)+k );
							coeff.Col( pascal(D+1,j)+pascal(i,j+1)+k ) = 0;
						}
				return coeff.Col(0);
			}

			Matrix<int> pascal_sym()
			{
				Matrix<int> pasc(D+2,ord+2);
				int i, j;

				for (i = 0; i <= D+1; ++i)
					for (j = 0; j <= ord+1; ++j)
						if (i == 0 || j == 0)
							pasc(i,j) = 0;
						else if (i == 1 || j == 1)
							pasc(i,j) = 1;
						else
							pasc(i,j) = pasc(i-1,j) + pasc(i,j-1);
				return pasc;
			}
	};

}

#endif
