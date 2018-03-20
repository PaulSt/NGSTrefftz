#include "helpers.hpp"

namespace ngfem
{

	// double Horner( Vector<double> a, int x )
	// {
	// 	int deg = a.Size()-1;
	// 	double result = a[deg];
	// 	for(int i=deg-1; i >= 0 ; --i)
	// 		result = result * x + a[i];
	// 	return result;
	// }
	// 
	// double MultiHorner( Matrix<int> multiind, Vector<int> klist, Vector<double> coeff, Vector<double> point)
	// {
	// 	int d = multiind.Width();
	// 	int M = multiind.Height()-1;
	// 	Vector<double> result(d+1);
	// 	result = 0;
	// 	int k;
	// 	result(0) = coeff(M);
	// 	double temp = 0;
	//
	// 	for(int n=1;n<=M;n++)
	// 	{
	// 		k = klist(M-n+1);
	// 		temp = 0;
	// 		for(int i=0;i<=k;i++)
	// 			temp += result(i);
	// 		result(k) = point(k-1)*temp;
	// 		result(0) = coeff(M-n);
	// 		for(int i=1;i<k;i++)
	// 			result(i) = 0;
	// 	}
	// 	temp=0;
	// 	for(int i=0;i<=d;i++)
	// 		temp += result(i);
	//
	// 	return temp;
	// }
	//
	// Vector<int> maxj(Matrix<int> multiind)
	// {
	// 	Vector<int> result(multiind.Height());
	// 	for(int i = multiind.Height()-1; i>=0; i--)
	// 	{
	// 		for(int d = multiind.Width()-1; d>=0; d--)
	// 		{
	// 			if(multiind(i,d)-multiind(i-1,d) != 0)
	// 			{
	// 				result(i) = d+1;
	// 				break;
	// 			}
	// 		}
	// 	}
	// 	return result;
	// }

	// user horner in calcshape:
	// Vector<int> klist(maxj(indices));
	// for(int b = 0; b < nbasis; b++) shape(b) = MultiHorner( indices, klist, basis.Row(b), cpoint);

}
