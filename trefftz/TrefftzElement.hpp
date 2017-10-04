#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include <l2hofefo.hpp>
#include"helpers.cpp"
#include "MappedElement.hpp"
//using namespace ngfem;

namespace ngfem
{
	template <int D, int ord>
	class TrefftzElement : public MappedElement
	{
		private:

			constexpr static int nbasis = BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1);

			constexpr static int npoly = BinCoeff(D + ord, ord);

			static const Mat<npoly, D, int> indices;

			Vec<D> elcenter;

			//static const Mat<nbasis, npoly,double> basis;
			static const Matrix<double> basis;

		protected:

			void static MakeIndices_inner(Mat<npoly, D, int> &indice, Vec<D, int> &numbers, int &count, int dim = D);

			constexpr static Mat<npoly, D, int> MakeIndices();

			constexpr static int IndexMap(Vec<D, int> index);

		public:

			TrefftzElement();

			TrefftzElement(Vec<D> aelcenter);

			virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

			virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;

			virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;

			using MappedElement::CalcShape;

			using MappedElement::CalcDShape;

			double ipow_ar(FlatVector<double> base, Vec<D, int> ex, float result = 1, int count = D) const;

			int GetNBasis() const;

			constexpr static Mat<nbasis, npoly,double> TrefftzBasis();

			TrefftzElement<D,ord> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}

	};
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement(py::module m);
#endif // NGS_PYTHON

#endif // FILE_TrefftzElement_HPP
