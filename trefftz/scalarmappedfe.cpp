#include "scalarmappedfe.hpp"
#include "h1lofe.hpp"
#include "l2hofe.hpp"


namespace ngfem
{

    void BaseScalarMappedElement ::
    CalcShape (const BaseMappedIntegrationPoint & mip,
               BareSliceVector<Complex> shape) const
    {
        CalcShape (mip, SliceVector<double> (ndof, 2*shape.Dist(), reinterpret_cast<double*> (&shape(0))));
        SliceVector<double> imag_part(ndof, 2*shape.Dist(), reinterpret_cast<double*> (&shape(0))+1);
        imag_part = 0.0;
    }

    void BaseScalarMappedElement ::
    CalcShape (const BaseMappedIntegrationRule & mir,
               SliceMatrix<> shape) const
    {
        for (int i = 0; i < mir.Size(); i++)
            CalcShape (mir[i], shape.Col(i));
    }

    void BaseScalarMappedElement ::
    CalcShape (const SIMD_BaseMappedIntegrationRule & mir,
               BareSliceMatrix<SIMD<double>> shape) const
    {
        cout<<"SIMD - CalcShape not overloaded"<< endl;
        throw ExceptionNOSIMD("SIMD - CalcShape not overloaded");
    }

    double BaseScalarMappedElement ::
    Evaluate (const BaseMappedIntegrationPoint & mip, BareSliceVector<double> x) const
    {
        VectorMem<20, double> shape(ndof);
        CalcShape (mip, shape);
        return InnerProduct (shape, x);
    }

    void BaseScalarMappedElement ::
    Evaluate (const BaseMappedIntegrationRule & mir, BareSliceVector<double> coefs, FlatVector<double> vals) const
    {
        for (size_t i = 0; i < mir.Size();i++) //.GetNIP(); i++)
            vals(i) = Evaluate (mir[i], coefs);
    }

    void BaseScalarMappedElement ::
    Evaluate (const SIMD_BaseMappedIntegrationRule & mir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const
    {
        cout<<"SIMD - Eval not overloaded"<< endl;
        throw ExceptionNOSIMD (string("Evaluate (simd) not implemented for class ")+typeid(*this).name());
    }

    void BaseScalarMappedElement ::
    Evaluate (const SIMD_BaseMappedIntegrationRule & mir, SliceMatrix<> coefs, BareSliceMatrix<SIMD<double>> values) const
    {
        for (size_t i = 0; i < coefs.Width(); i++)
            Evaluate (mir, coefs.Col(i), values.Row(i));
    }

    void BaseScalarMappedElement ::
    Evaluate (const BaseMappedIntegrationRule & mir, SliceMatrix<> coefs, SliceMatrix<> values) const
    {
        VectorMem<100> shapes(coefs.Height());
        for (size_t i = 0; i < mir.Size(); i++)
        {
            CalcShape (mir[i], shapes);
            values.Row(i) = Trans(coefs) * shapes;
        }
    }

    void BaseScalarMappedElement ::
    EvaluateTrans (const BaseMappedIntegrationRule & mir, FlatVector<double> vals, BareSliceVector<double> coefs) const
    {
        VectorMem<20, double> shape(ndof);
        coefs.Range(0,ndof) = 0.0;
        for (int i = 0; i < mir.Size(); i++) //GetNIP()
        {
            CalcShape (mir[i], shape);
            coefs.Range(0,ndof) += vals(i) * shape;
        }
    }

    void BaseScalarMappedElement ::
    AddTrans (const SIMD_BaseMappedIntegrationRule & mir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const
    {
        cout<<"SIMD - AddTrans not overloaded"<< endl;
        throw ExceptionNOSIMD (string("AddTrans (simd) not implemented for class ")+typeid(*this).name());
    }

    void BaseScalarMappedElement ::
    AddTrans (const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const
    {
        for (int i = 0; i < coefs.Width(); i++)
            AddTrans (mir, values.Row(i), coefs.Col(i));
    }

    void BaseScalarMappedElement ::
    CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir,
                      BareSliceMatrix<SIMD<double>> dshapes) const
    {
        cout<<"SIMD - CalcMappedDShape not overloaded"<< endl;
        throw ExceptionNOSIMD("SIMD - CalcDShape not overloaded");
    }

    void BaseScalarMappedElement ::
    EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    {
        cout<<"SIMD - EvaluateGrad not overloaded"<< endl;
        throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for class ")+typeid(*this).name());
    }

    //void BaseScalarMappedElement ::
    //EvaluateGrad (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
    //{
    //throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for class ")+typeid(*this).name());
    //}

    void BaseScalarMappedElement ::
    AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                  BareSliceVector<> coefs) const
    {
        cout<<"SIMD - AddTransGrad not overloaded"<< endl;
        throw ExceptionNOSIMD (string("AddGradTrans (simd) not implemented for class ")+typeid(*this).name());
    }



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    template<int D>
    string ScalarMappedElement<D> :: ClassName() const
    {
        return "ScalarMappedElement";
    }


    template<int D>
    void ScalarMappedElement<D> :: CalcDShape (const BaseMappedIntegrationRule & mir,
                                               BareSliceMatrix<> dshapes) const
    {
        for (int i = 0; i < mir.Size(); i++)
            CalcDShape (mir[i], dshapes.Cols(i*D,(i+1)*D));
    }


    template<int D>
    void ScalarMappedElement<D> ::
    CalcMappedDShape (const BaseMappedIntegrationPoint & mip,
                      BareSliceMatrix<> dshape) const
    {
        //no mapping - no inner derivative
        CalcDShape (mip, dshape);
    }

    template<int D>
    void ScalarMappedElement<D> ::
    CalcMappedDShape (const BaseMappedIntegrationRule & mir,
                      SliceMatrix<> dshapes) const
    {
        for (int i = 0; i < mir.Size(); i++)
            CalcMappedDShape (mir[i], dshapes.Cols(i*D,(i+1)*D));
    }


    template<int D>
    Vec<D> ScalarMappedElement<D> ::
    EvaluateGrad (const BaseMappedIntegrationPoint & ip, BareSliceVector<double> x) const
    {
        MatrixFixWidth<D> dshape(ndof);
        CalcDShape (ip, dshape);
        Vec<D> grad = Trans (dshape) * x;
        return grad;
    }


    template<int D>
    void ScalarMappedElement<D> ::
    EvaluateGrad (const BaseMappedIntegrationRule & ir, BareSliceVector<double> coefs, FlatMatrixFixWidth<D,double> vals) const
    {
        for (size_t i = 0; i < ir.Size(); i++)
            vals.Row(i) = EvaluateGrad (ir[i], coefs);
    }



    template<int D>
    void ScalarMappedElement<D> ::
    EvaluateGradTrans (const BaseMappedIntegrationRule & ir, FlatMatrixFixWidth<D,double> vals,
                       BareSliceVector<double> coefs) const
    {
        MatrixFixWidth<D> dshape(ndof);
        coefs.Range(0,ndof) = 0.0;
        for (int i = 0; i < ir.Size(); i++)
        {
            CalcDShape (ir[i], dshape);
            coefs.Range(0,ndof) += dshape * vals.Row(i);
        }
    }


    template<int D>
    void ScalarMappedElement<D> ::
    EvaluateGradTrans (const BaseMappedIntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const
    {
#ifndef __CUDA_ARCH__
        cout << "EvalGradTrans not overloaded" << endl;
#endif
    }

    //template<int D>
    //void ScalarMappedElement<D> ::
    //GetPolOrders (FlatArray<PolOrder<D> > orders) const
    //{
//#ifndef __CUDA_ARCH__
        //throw Exception (string ("GetPolOrders not implemnted for element") + ClassName());
//#endif
    //}


    template<int D>
    void ScalarMappedElement<D> ::
    CalcDShape (const SIMD_BaseMappedIntegrationRule & smir, BareSliceMatrix<SIMD<double>> dshape) const
    {
        cout<<"SIMD - CalcDShape not overloaded"<< endl;
        throw ExceptionNOSIMD("SIMD - CalcDShape not overloaded");
    }



    template<>
    void ScalarMappedElement<4> :: CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip,
                                                BareSliceMatrix<> hddshape) const
    {;}

    template<int D>
    void ScalarMappedElement<D> :: CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip,
                                                BareSliceMatrix<> hddshape) const
    {
        auto & mip = static_cast<const MappedIntegrationPoint<D,D> &> (bmip);
        int nd = GetNDof();
        auto ddshape = hddshape.AddSize(nd, D*D);
        double eps = 1e-7;
        MatrixFixWidth<(D)> dshape1(nd), dshape2(nd);
        const ElementTransformation & eltrans = mip.GetTransformation();

        for (int i = 0; i < D; i++)
        {
            IntegrationPoint ip1 = mip.IP();
            IntegrationPoint ip2 = mip.IP();
            ip1(i) -= eps;
            ip2(i) += eps;
            MappedIntegrationPoint<D,D> mip1(ip1, eltrans);
            MappedIntegrationPoint<D,D> mip2(ip2, eltrans);

            //cout << bmip.GetPoint() << endl;
            //bmip.GetPoint()(i) -= eps;
            //cout << bmip.GetPoint() << endl;
            CalcMappedDShape (mip1, dshape1);
            //bmip.GetPoint()(i) += 2*eps;
            CalcMappedDShape (mip2, dshape2);
            //bmip.GetPoint()(i) -= eps;

            ddshape.Cols(D*i,D*(i+1)) = (0.5/eps) * (dshape2-dshape1);
        }

        for (int j = 0; j < D; j++)
        {
            for (int k = 0; k < nd; k++)
                for (int l = 0; l < D; l++)
                    dshape1(k,l) = ddshape(k, l*D+j);

            dshape2 = dshape1 * mip.GetJacobianInverse();

            for (int k = 0; k < nd; k++)
                for (int l = 0; l < D; l++)
                    ddshape(k, l*D+j) = dshape2(k,l);
        }
    }




    template class ScalarMappedElement<1>;
    template class ScalarMappedElement<2>;
    template class ScalarMappedElement<3>;
    template class ScalarMappedElement<4>;

}
