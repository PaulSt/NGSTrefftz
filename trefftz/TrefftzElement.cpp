#include <fem.hpp>
#include "TrefftzElement.hpp"
//#include "helpers.cpp"

namespace ngfem
{

	template<int D, int ord>
	const Mat<TrefftzElement<D,ord>::npoly, D, int> TrefftzElement<D,ord> :: indices = MakeIndices();

	//template<int D, int ord>
	//const Mat<TrefftzElement<D,ord>::nbasis, TrefftzElement<D,ord>::npoly,double> TrefftzElement<D,ord> :: basis = TrefftzBasis();


	template <int D, int ord>
  void TrefftzElement<D,ord> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                  BareSliceVector<> shape) const
  {
		Vec<npoly,float> polynomial;

		for(int i=0; i<npoly; i++)//loop over indices
		{
			polynomial(i) = ipow_ar(mip.GetPoint(),indices.Row(i));
		}

		Vec<nbasis,double> tempshape = basis * polynomial;
		for(int i = 0; i< nbasis; i++) shape(i) = tempshape(i);

/*
		FlatVector<double> point = mip.GetPoint();
		shape(0) = point(0) * point(1);
*/
	}

	template <int D, int ord>
	void TrefftzElement<D,ord> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
																		SliceMatrix<> dshape) const
	{
		/*
		array<int, D> tempexp;
		int coeff;
		for(int l=0;l<nbasis;l++) //loop over basis functions
		{
			for(int d=0;d<D;d++)  //loop over derivatives/dimensions
			{
				for(int i=0;i<BinCoeff(D + ord, ord);i++)//loop over indices
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


	template <int D, int ord>
	double TrefftzElement<D,ord> :: ipow_ar(FlatVector<double> base, Vec<D, int> ex, float result, int count) const
	{
		return count == 0 ? result : ipow_ar( base, ex, pow(base(count-1),ex(count-1)) * result, count-1 );
	}

	template <int D, int ord>
	int TrefftzElement<D,ord> :: GetNBasis() const
	{
		return nbasis;
	}

}


#ifdef NGS_PYTHON
void ExportTrefftzElement(py::module m)
{
  using namespace ngfem;
  py::class_<TrefftzElement<3,3>, shared_ptr<TrefftzElement<3,3>>, FiniteElement>
    (m, "TrefftzElement", "Trefftz space for wave eq")
    //.def(py::init<>())
		.def(py::init<>())
		.def("TrefftzBasis", &TrefftzElement<3,3>::TrefftzBasis)
		//.def("CalcShape", &TrefftzElement<2>::CalcShape)
		.def("GetNBasis", &TrefftzElement<3,3>::GetNBasis)
		;
}
#endif // NGS_PYTHON
