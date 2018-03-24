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


	class HornerScheme
	{
		private:
			int D;
			int ord;
			int npoly;
			Matrix<int> pascal;
		public:
			HornerScheme(int aD, int aord): D(aD), ord(aord), npoly(BinCoeff(D + ord, ord)), pascal(pascal_sym())
			{;}

			void MultiHorner(Vector<double> coeff, Vector<double> point) const
			{
				for(int j=ord;j>0;j--)
					for(int i=1;i<=D;i++)
					 	for(int k=1;k<=pascal(j-1,i-1);k++){
if(coeff.Size() < pascal(j-1,D)+pascal(j,i-2)+k-1){
							// cout << pascal << endl;
							// cout << coeff.Size() << endl;

							// cout << "D: " << D << " ord: " << ord <<" j: " << j << " i: " << i << " k: " << k<< " index1: " <<pascal(j-2,D)+k << " index2: " <<  pascal(j-1,D)+pascal(j,i-2)+k << endl<<  endl;
						}// coeff( pascal(j-2,D)+k ) += point( D-i+1 ) * coeff( pascal(j-1,D)+pascal(j,i-2)+k );
							// coeff( pascal(j-1,D)+pascal(j,i-2)+k ) = 0;
							// cout << coeff <<  endl;
						}
			}

			Matrix<int> pascal_sym() {
				Matrix<int> pasc(ord+1,D+1);
				int i, j;

				for (i = 0; i <= ord; ++i)
					for (j = 0; j <= D; ++j)
						if (i == 0 || j == 0)
							pasc(i,j) = 1;
						else
							pasc(i,j) = pasc(i-1,j) + pasc(i,j-1);
				return pasc;
			}
	};

}

#endif
