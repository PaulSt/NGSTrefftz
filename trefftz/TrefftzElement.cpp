#include <fem.hpp>
#include "TrefftzElement.hpp"
//#include "helpers.cpp"

namespace ngfem
{

	template <int D, int order>
  void TrefftzElement<D,order> :: CalcShape (const IntegrationPoint & ip,
                                  BareSliceVector<> shape) const
  {

	}

	template <int D, int order>
	void TrefftzElement<D,order> :: CalcDShape (const IntegrationPoint & ip,
																		SliceMatrix<> dshape) const
	{

	}

	template <int D, int order>
  void TrefftzElement<D,order> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                  Vector<double> &shape) const
  {
		Vector<float> polynomial( BinCoeff(D+1 + order, order) );

		for(int i=0;i<BinCoeff(D+1 + order, order);i++)//loop over indices
		{
			polynomial(i) = ipow_ar(mip.GetPoint(),indices[i]);
		}
		shape = basis * polynomial;

		//FlatVector<double> point = mip.GetPoint();
		//shape(0) = point(0) * point(1);
	}

	template <int D, int order>
	void TrefftzElement<D,order> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
																		SliceMatrix<> dshape) const
	{
		/*
		array<int, D+1> tempexp;
		int coeff;
		for(int l=0;l<nbasis;l++) //loop over basis functions
		{
			for(int d=0;d<D;d++)  //loop over derivatives/dimensions
			{
				for(int i=0;i<BinCoeff(D+1 + order, order);i++)//loop over indices
				{
					if(indices[i][d+1] == 0) continue;
					else
					{
						tempexp = indices[i];
						dshape(l,d) += tempexp[d+1]-- * basisFunctions[l].get(indices[i]) * ipow_ar(mip.GetPoint(),tempexp,1,D+1);
					}
				}
			}
		}
		*/
	}


	template <int D, int order>
	void TrefftzElement<D,order> :: TrefftzBasis()
	{
		basis = 0;
		for(int l=0;l<nbasis;l++) //loop over basis functions
		{
			for(int i=0;i<indices.size();i++)//loop over indices BinCoeff(D+1 + order, order)
			{
				int k = indices[i][0];
				if(k > 1 )
				{
					for(int m=1;m<=D;m++) //rekursive sum
					{
						array<int, D+1> get_coeff = indices[i];
						get_coeff[0] = get_coeff[0] - 2;
						get_coeff[m] = get_coeff[m] + 2;
						basis( l, IndexMap(indices[i]) ) += (indices[i][m]+1) * (indices[i][m]+2) * basis(l,IndexMap(get_coeff) );
					}
					basis( l, IndexMap(indices[i]) ) *= 1.0/(k * (k-1));
				}
				else if(k == 0 ) //time=0
				{
					basis( l, IndexMap(indices[l]) ) = 1.0;
					i += BinCoeff(D + order, order) + BinCoeff(D + order-1, order-1);
				}
			}
		}
		//cout << "basis: \n" << basis << endl;
	}


	template <int D, int order>
	void TrefftzElement<D,order> :: MakeIndices(array<int, D+1> &numbers, int &count, int dim)
	{
		if (dim>0)
		{
			for(int i=0;i<=order;i++)
			{
				numbers[numbers.size() - dim]=i;
				MakeIndices(numbers,count,dim-1) ;
			}
		}
		else
		{
			int sum=0;
			for(int i=0;i<numbers.size();i++)
			{
				sum += numbers[i];
			}
			if(sum<=order){
				indices[count++] = numbers;
				/*
				cout << IndexMap(indices[count-1]) << ": ";
				for(int i=0;i<numbers.size();i++)
				{
					cout << indices[count-1][i] <<" ";
				}
				cout << "\n";
				*/
			}
		}
	}


	template <int D, int order>
	float TrefftzElement<D,order> :: ipow_ar(IntegrationPoint base, array<int, D+1> ex, float result, int count) const
	{
		return count == 0 ? result * pow(base(count),ex[count]) : ipow_ar( base, ex, result * pow(base(count),ex[count]), count-1 );
	}

	template <int D, int order>
	double TrefftzElement<D,order> :: ipow_ar(FlatVector<double> base, array<int, D+1> ex, float result, int count) const
	{
		//cout << "run: " << count << " result: "<< result << " calc: " << base(count-1) << " hoch " << ex[count-1] << " berechnung: " << pow(base(count-1),ex[count-1]) * result << endl;
		return count == 0 ? result : ipow_ar( base, ex, pow(base(count-1),ex[count-1]) * result, count-1 );
	}


	template <int D, int order>
	int TrefftzElement<D,order> :: GetNBasis() const
	{
		return nbasis;
	}

	template <int D, int order>
	inline int TrefftzElement<D,order> :: IndexMap(array<int, D+1> index)
	{
		int sum=0;
		int temp_size = 0;
		for(int d=0;d<D+1;d++){
			for(int p=0;p<index[d];p++){
				sum+=BinCoeff(D - d + order - p - temp_size, order - p - temp_size);
			}
			temp_size+=index[d];
		}
		return sum;
	}

}


#ifdef NGS_PYTHON
void ExportTrefftzElement(py::module m)
{
  using namespace ngfem;
  py::class_<TrefftzElement<2,3>, shared_ptr<TrefftzElement<2,3>>, BaseScalarFiniteElement>
    (m, "TrefftzElement", "Trefftz space for wave eq")
    //.def(py::init<>())
		.def(py::init<>())
		.def("TrefftzBasis", &TrefftzElement<2,3>::TrefftzBasis)
		//.def("CalcShape", &TrefftzElement<2>::CalcShape)
		.def("GetNBasis", &TrefftzElement<2,3>::GetNBasis)
		;
}
#endif // NGS_PYTHON
