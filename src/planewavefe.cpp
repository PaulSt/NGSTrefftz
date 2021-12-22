#include "planewavefe.hpp"

namespace ngfem
{

    template<>
    Vec<2> PlaneWaveElement<2> :: GetDirection(int i) const
    {
        return { cos(2.0*M_PI*i/this->ndof), sin(2.0*M_PI*i/this->ndof) };
    }

    template<>
    void PlaneWaveElement<2> :: CalcShape (const BaseMappedIntegrationPoint & mip,
                                        BareSliceVector<Complex> shape) const
    {
        Vec<2> cpoint = mip.GetPoint();
        cpoint -= elcenter;

        for (int i=0; i<this->ndof; ++i)
        {
            Vec<2> dir = GetDirection(i);
            shape(i) = 0.0;
            shape(i) = exp(1i*(cpoint[0]*dir[0] + cpoint[1]*dir[1])*conj);
        }
    }

    template<>
    void PlaneWaveElement<2> :: CalcDShape (const BaseMappedIntegrationPoint & mip, BareSliceMatrix<Complex> dshape) const
    {
        Vec<2> cpoint = mip.GetPoint();
        cpoint -= elcenter;

        for(int d=0;d<2;d++)
            for (int i=0; i<this->ndof; ++i)
            {
                Vec<2> dir = GetDirection(i);
                dshape(i,d) = 0.0;
                dshape(i,d) = 1i*dir[d]*conj*exp(1i*(cpoint[0]*dir[0] + cpoint[1]*dir[1])*conj);
            }
    }

    template<int D>
    Complex PlaneWaveElement<D> ::
    EvaluateComplex (const BaseMappedIntegrationPoint & mip, BareSliceVector<Complex> x) const
    {
        VectorMem<20, Complex> shape(this->ndof);
        CalcShape (mip, shape);
        return InnerProduct (shape, x);
    }

    template<int D>
    void PlaneWaveElement<D> ::
    Evaluate (const BaseMappedIntegrationRule & mir, BareSliceVector<Complex> coefs, FlatVector<Complex> vals) const
    {
        for (size_t i = 0; i < mir.Size();i++) //.GetNIP(); i++)
            vals(i) = EvaluateComplex (mir[i], coefs);
    }

    template<int D>
    Vec<D,Complex> PlaneWaveElement<D> ::
    EvaluateGradComplex (const BaseMappedIntegrationPoint & ip, BareSliceVector<Complex> x) const
    {
        MatrixFixWidth<D,Complex> dshape(this->ndof);
        CalcDShape (ip, dshape);
        Vec<D,Complex> grad = Trans (dshape) * x;
        return grad;
    }


    template class PlaneWaveElement<2>;
}
