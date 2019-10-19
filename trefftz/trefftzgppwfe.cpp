#include "trefftzgppwfe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"
#include "helpers.hpp"

#include <ctime>

namespace ngfem
{
    template<int D>
    TrefftzGppwFE<D> :: TrefftzGppwFE(int aord, float ac, Vec<D> aelcenter, double aelsize, ELEMENT_TYPE aeltype, int abasistype)
    : ScalarMappedElement<D+1>(2*aord+1, aord),
    ord(aord),
    c(ac),
    npoly(BinCoeff(D + ord, ord)),
    elcenter(aelcenter),
    elsize(aelsize),
    eltype(aeltype),
    basistype(abasistype)
    {;}


    template<>
    void TrefftzGppwFE<1> :: CalcShape (const SIMD_BaseMappedIntegrationRule & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;
            // calc 1 dimensional monomial basis
            STACK_ARRAY(SIMD<double>, mem, 2*(ord+1)+1);
            for(size_t d=0;d<2;d++)
            {
                SIMD<double> evalpoint = pow(-1,d)*cpoint[0]-cpoint[1];
                Monomial (ord, evalpoint, &mem[d*(ord+1)]);
            }
            for (int i=0; i<this->ndof; ++i)
            {
                shape(i,imip) = i<=ord ? mem[i] : mem[i+1];
            }
        }
    }

    template<>
    void TrefftzGppwFE<2> :: CalcShape (const SIMD_BaseMappedIntegrationRule & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {cout << "dim not implemented" << endl;}

    template<>
    void TrefftzGppwFE<3> :: CalcShape (const SIMD_BaseMappedIntegrationRule & smir,
                                        BareSliceMatrix<SIMD<double>> shape) const
    {
    }



    template<>
    void TrefftzGppwFE<1> :: CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
                                         BareSliceMatrix<SIMD<double>> dshape) const
    {
        for (int imip = 0; imip < smir.Size(); imip++)
        {
            Vec<2,SIMD<double>> cpoint = smir[imip].GetPoint();
            cpoint -= elcenter; cpoint *= (2.0/elsize); cpoint[1] *= c;
            // calc 1 dimensional monomial basis
            STACK_ARRAY(SIMD<double>, mem, 2*(ord+1)+1+2);
            for(size_t d=0;d<2;d++)
            {
                SIMD<double> evalpoint = pow(-1,d)*cpoint[0]-cpoint[1];
                Monomial (ord, evalpoint, &mem[d*(ord+1)]);
            }
            for(int d=0;d<2;d++)
            {
                for (int i=0; i<this->ndof; ++i)
                {
                    dshape(i*2+d,imip) = pow(-1,(i<=ord && d==0)) * (i<=ord ? i : i-ord)  * (i<=ord ? mem[i-1] : mem[i])  * (d==1 ? c : 1) * (2.0/elsize);
                }
            }
        }
    }

    template<>
    void TrefftzGppwFE<2> :: CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
                                         BareSliceMatrix<SIMD<double>> dshape) const
    {}

    template<>
    void TrefftzGppwFE<3> :: CalcDShape (const SIMD_BaseMappedIntegrationRule & smir,
                                         BareSliceMatrix<SIMD<double>> dshape) const
    {
    }


    /////////////// non-simd


    template<>
    void TrefftzGppwFE<1> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {cout << "dim not implemented" << endl;}

    template<>
    void TrefftzGppwFE<2> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
    }

    template<>
    void TrefftzGppwFE<3> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<> shape) const
    {
    }



    template<>
    void TrefftzGppwFE<1> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {cout << "dim not implemented" << endl;}

    template<>
    void TrefftzGppwFE<2> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {
    }

    template<>
    void TrefftzGppwFE<3> :: CalcDShape (const BaseMappedIntegrationPoint & mip,
                                         SliceMatrix<> dshape) const
    {
    }


    template class TrefftzGppwFE<1>;
    template class TrefftzGppwFE<2>;
    template class TrefftzGppwFE<3>;


    template<int D>
    const CSR* TrefftzGppwBasis<D> :: TB(int ord)
    {
        const CSR* tb =& tbstore[ord];
        return tb;
    }

    template<int D>
    void TrefftzGppwBasis<D> :: CreateTB(int ord, int basistype)
    {
        cout << "creating tp store for order " << ord << endl;
    }


    template<int D>
    void TrefftzGppwBasis<D> :: TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D, int> coeffnum, int basis, int dim, int &tracker, int basistype)
    {
    }

    template<int D>
    int TrefftzGppwBasis<D> :: IndexMap2(Vec<D, int> index, int ord)
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

    template class TrefftzGppwBasis<1>;
    template class TrefftzGppwBasis<2>;
    template class TrefftzGppwBasis<3>;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//#ifdef NGS_PYTHON
//void ExportTrefftzElement(py::module m)
//{
//// py::class_<TrefftzGppwFE<3>, shared_ptr<TrefftzGppwFE<3>>, FiniteElement>
//// 	(m, "TrefftzWaveFE3", "Trefftz space for wave eq")
//// 	.def(py::init<>())
//// 	;
//// py::class_<TrefftzGppwFE<2>, shared_ptr<TrefftzGppwFE<2>>, FiniteElement>
//// 	(m, "TrefftzWaveFE2", "Trefftz space for wave eq")
//// 	.def(py::init<>())
//// 	;
//// py::class_<TrefftzGppwFE<1>, shared_ptr<TrefftzGppwFE<1>>, FiniteElement>
//// 	(m, "TrefftzWaveFE1", "Trefftz space for wave eq")
//// 	.def(py::init<>())
//// 	;
//}
//#endif // NGS_PYTHON
