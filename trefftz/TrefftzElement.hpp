#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include "helpers.cpp"
#include "ScalarMappedElement.hpp"

namespace ngfem
{
	template <int D>
	class T_TrefftzElement : public ScalarMappedElement<D>
	{
		private:
			const int ord;
			const int nbasis;
			const int npoly;
			const Matrix<int> indices;
			const Matrix<double> basis;
			Vec<D> elcenter=0; float elsize=1; float c=1;
			ELEMENT_TYPE eltype;

		public:
			T_TrefftzElement();
			T_TrefftzElement(int aord);
			T_TrefftzElement(int aord, ELEMENT_TYPE aeltype);

			virtual ELEMENT_TYPE ElementType() const { return eltype; }

			using ScalarMappedElement<D>::CalcShape;
			virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;

			using ScalarMappedElement<D>::CalcDShape;
			virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;

			int GetNBasis() const { return nbasis; }

			T_TrefftzElement<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
			T_TrefftzElement<D> * SetElSize(float aelsize) {elsize = aelsize; return this;}
			T_TrefftzElement<D> * SetWavespeed(float ac) {c = ac; return this;}

		protected:
			FlatVector<double> ShiftPoint(FlatVector<double> point) const;

			constexpr void MakeIndices_inner(Matrix<int> &indice, Vec<D, int> &numbers, int &count, int dim = D);
			constexpr Matrix<int> MakeIndices();

			constexpr int IndexMap(Vec<D, int> index) const;
			double ipow_ar(FlatVector<double> base, Vec<D, int> ex, double result = 1, int count = D-1) const;
			double ipowD_ar(int der, FlatVector<double> base, Vec<D, int> ex, double result = 1, int count = D-1) const;
			constexpr Matrix<double> TrefftzBasis() const;
	};
}


#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement(py::module m);
#endif // NGS_PYTHON


#endif // FILE_TrefftzElement_HPP
