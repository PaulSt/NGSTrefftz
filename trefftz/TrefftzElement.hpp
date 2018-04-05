#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include "helpers.hpp"
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
			const int basistype;
			Vec<D> elcenter=0;
			double elsize=1;
			float c = 1.0;
			ELEMENT_TYPE eltype;
			HornerScheme<D> horner;

		public:
			// T_TrefftzElement();
			T_TrefftzElement(int aord = 1, float ac = 1.0, ELEMENT_TYPE aeltype = ET_TRIG, int abasistype = 0);

			virtual ELEMENT_TYPE ElementType() const { return eltype; }

			using ScalarMappedElement<D>::CalcShape;
			virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;

			using ScalarMappedElement<D>::CalcDShape;
			virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;

			int GetNBasis() const { return nbasis; }

			T_TrefftzElement<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
			T_TrefftzElement<D> * SetElSize(double aelsize) {elsize = aelsize; return this;}
			// T_TrefftzElement<D> * SetWavespeed(float ac) {c = ac; return this;}

		protected:
			Vector<double> ShiftPoint(Vector<double> point) const;

			constexpr void MakeIndices_inner(Matrix<int> &indice, Vec<D, int> &numbers, int &count, int ordr, int dim = D);
			Matrix<int> GetIndices();

			constexpr int IndexMap(Vec<D, int> index) const;
			Matrix<double> TrefftzBasis() const;
			Matrix<double> GetDerTrefftzBasis(int der) const;
	};
}


#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement(py::module m);
#endif // NGS_PYTHON


#endif // FILE_TrefftzElement_HPP
