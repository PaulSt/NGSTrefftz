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
			ELEMENT_TYPE eltype;

		protected:
			void static MakeIndices_inner(Mat<npoly, D, int> &indice, Vec<D, int> &numbers, int &count, int dim = D);
			constexpr static Mat<npoly, D, int> MakeIndices();
			constexpr static int IndexMap(Vec<D, int> index);

		public:
			TrefftzElement()
				: ScalarMappedElement<D>(nbasis,ord), //BaseScalarMappedElement(nbasis,ord)
					eltype(ET_TRIG) {;}
			TrefftzElement(ELEMENT_TYPE aeltype)
				: ScalarMappedElement<D>(nbasis,ord),
					eltype(aeltype) {;}

			virtual ELEMENT_TYPE ElementType() const { return eltype; }

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template <int D>
	class T_TrefftzElement : public ScalarMappedElement<D>
	{
		private:
			const int ord;
			const int nbasis = BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1);
			const int npoly = BinCoeff(D + ord, ord);
			const Matrix<int> indices = MakeIndices();
			const Matrix<double> basis = TrefftzBasis();
			Vec<D> elcenter=0; float elsize=1; float c=1;
			ELEMENT_TYPE eltype;

		public:
			T_TrefftzElement()
				: ord(1),
					ScalarMappedElement<D>(3, 1),
					eltype(ET_TRIG){;}
			T_TrefftzElement(int aord)
				:	ord(aord),
					ScalarMappedElement<D>(BinCoeff(D-1 + aord, aord) + BinCoeff(D-1 + aord-1, aord-1), aord),
					eltype(ET_TRIG) {;}
			T_TrefftzElement(int aord, ELEMENT_TYPE aeltype)
				: ord(aord),
					ScalarMappedElement<D>(BinCoeff(D-1 + aord, aord) + BinCoeff(D-1 + aord-1, aord-1), aord),
					eltype(aeltype) {;}

			virtual ELEMENT_TYPE ElementType() const { return eltype; }

			using ScalarMappedElement<D>::CalcShape;
			virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;

			using ScalarMappedElement<D>::CalcDShape;
			virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, SliceMatrix<> dshape) const;

			int GetNBasis() const { return nbasis; }

			T_TrefftzElement<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
			T_TrefftzElement<D> * SetElSize(float aelsize) {elsize = aelsize; return this;}
			T_TrefftzElement<D> * SetWavespeed(float ac) {c = ac; return this;}

			FlatVector<double> ShiftPoint(FlatVector<double> point) const
			{ point -= elcenter; point *= (1.0/elsize); point(0) *= c; return point; }

			double ipow_ar(FlatVector<double> base, Vec<D, int> ex, double result = 1, int count = D-1) const;

			double ipowD_ar(int der, FlatVector<double> base, Vec<D, int> ex, double result = 1, int count = D-1) const;

			constexpr Matrix<double> TrefftzBasis() const
			{
				Matrix<double> temp_basis(nbasis,npoly);
				temp_basis = 0;
				for(int l=0;l<nbasis;l++) //loop over basis functions
				{
					for(int i=0;i<npoly;i++)//loop over indices BinCoeff(D + ord, ord)
					{
						int k = indices(i,0);
						if(k > 1 )
						{
							for(int m=1;m<=D-1;m++) //rekursive sum
							{
								Vec<D, int> get_coeff = indices.Row(i);
								get_coeff[0] = get_coeff[0] - 2;
								get_coeff[m] = get_coeff[m] + 2;
								temp_basis( l, IndexMap(indices.Row(i)) ) += (indices(i,m)+1) * (indices(i,m)+2) * temp_basis(l,IndexMap(get_coeff) );
							}
							temp_basis( l, IndexMap(indices.Row(i)) ) *= 1.0/(k * (k-1));
						}
						else if(k == 0 ) //time=0
						{
							temp_basis( l, IndexMap(indices.Row(l)) ) = 1.0;
							i += nbasis-1;
						}
					}
				}
				return temp_basis;
			}

		protected:
			constexpr void MakeIndices_inner(Matrix<int> &indice, Vec<D, int> &numbers, int &count, int dim = D)
			{
				if (dim>0)
				{
					for(int i=0;i<=ord;i++)
					{
						numbers(D - dim)=i;
						MakeIndices_inner(indice,numbers,count,dim-1) ;
					}
				}
				else
				{
					int sum=0;
					for(int i=0;i<D;i++)
					{
						sum += numbers(i);
					}
					if(sum<=ord){
						indice.Row(count++) = numbers;
					}
				}
			}

			constexpr Matrix<int> MakeIndices()
			{
				Matrix<int> indice(npoly,D);
				Vec<D, int>  numbers = 0;
				int count = 0;
				MakeIndices_inner(indice, numbers, count);
				return indice;
			}

			constexpr int IndexMap(Vec<D, int> index) const
			{
				int sum=0;
				int temp_size = 0;
				for(int d=0;d<D;d++){
					for(int p=0;p<index(d);p++){
						sum+=BinCoeff(D-1 - d + ord - p - temp_size, ord - p - temp_size);
					}
					temp_size+=index(d);
				}
				return sum;
			}
	};


}




#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement(py::module m);
#endif // NGS_PYTHON


#endif // FILE_TrefftzElement_HPP
