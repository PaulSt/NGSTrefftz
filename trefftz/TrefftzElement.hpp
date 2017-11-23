#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include "helpers.cpp"
#include "ScalarMappedElement.hpp"

namespace ngfem
{
	template <int D, int ord>
	class TrefftzElement : public ScalarMappedElement<D>
	{
		private:
			constexpr static int nbasis = BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1);
			constexpr static int npoly = BinCoeff(D + ord, ord);
			static const Mat<npoly, D, int> indices;
			//static const Mat<nbasis, npoly,double> basis;
			static const Matrix<double> basis;
			Vec<D> elcenter=0; float elsize=1; float c=1;

		protected:
			void static MakeIndices_inner(Mat<npoly, D, int> &indice, Vec<D, int> &numbers, int &count, int dim = D);
			constexpr static Mat<npoly, D, int> MakeIndices();
			constexpr static int IndexMap(Vec<D, int> index);

		public:
			TrefftzElement(): ScalarMappedElement<D>(nbasis,ord) { ;} //BaseScalarMappedElement(nbasis,ord) { ;	}//

			virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

			using ScalarMappedElement<D>::CalcShape;
			virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;

			using ScalarMappedElement<D>::CalcDShape;
			virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;

			constexpr static Mat<nbasis, npoly,double> TrefftzBasis();

			int GetNBasis() const { return nbasis; }

			TrefftzElement<D,ord> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
			TrefftzElement<D,ord> * SetElSize(float aelsize) {elsize = aelsize; return this;}
			TrefftzElement<D,ord> * SetWavespeed(float ac) {c = ac; return this;}

			FlatVector<double> ShiftPoint(FlatVector<double> point) const
			{ point -= elcenter; point *= (1.0/elsize); point(0) *= c; return point; }

			double ipow_ar(FlatVector<double> base, Vec<D, int> ex, double result = 1, int count = D-1) const;

			double ipowD_ar(int der, FlatVector<double> base, Vec<D, int> ex, double result = 1, int count = D-1) const;
	};



//////////////////////////////////////////////////////////////////////////////////////////////////////////


	template <int D, int ord>
	class TrefftzPWElement : public ScalarMappedElement<D>
	{
		private:
			constexpr static int nbasis = 2*ord+1;
			// constexpr static int npoly = ;
			Vec<D> elcenter=0; float c=1;
			static const Matrix<double> directions;

		public:
			TrefftzPWElement(): ScalarMappedElement<D>(nbasis,ord) { ; } //BaseScalarMappedElement(nbasis,ord) { ;	}//

			virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

			using ScalarMappedElement<D>::CalcShape;
			virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;

			using ScalarMappedElement<D>::CalcDShape;
			virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;

			int GetNBasis() const { return nbasis; }

			constexpr static Mat<nbasis, D,double> MakeDirections();

			TrefftzPWElement<D,ord> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
			TrefftzPWElement<D,ord> * SetWavespeed(float ac) {c = ac; return this;}

	};
}




#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement(py::module m);
#endif // NGS_PYTHON


#endif // FILE_TrefftzElement_HPP
