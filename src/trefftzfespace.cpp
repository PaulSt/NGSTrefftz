#include <comp.hpp>    // provides FESpace, ...
#include <python_comp.hpp>

#include "trefftzfespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{

    TrefftzFESpace :: TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
        : FESpace (ama, flags)
    {
        type="trefftzfespace";

        D = ma->GetDimension();

        this->dgjumps = true;
        order = int(flags.GetNumFlag ("order", 1));
        c = flags.GetNumFlag ("wavespeed", 1);
        wavespeedcf = make_shared<ConstantCoefficientFunction>(c);
        basistype = flags.GetNumFlag ("basistype", 0);
        useshift = flags.GetNumFlag("useshift",1);
        usescale = flags.GetNumFlag("usescale",1);
        DefineNumListFlag("eq");
        eqtyp = flags.GetStringFlag("eq");

        if(eqtyp=="fowave" || eqtyp=="foqtwave")
            local_ndof = (D)*BinCoeff(D-1 + order, D-1);
        else if(eqtyp=="wave")
            local_ndof = BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1) - (eqtyp=="fowave_reduced");
        else if(eqtyp=="laplace")
            local_ndof = BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1);
        else if(eqtyp=="helmholtz" || eqtyp=="helmholtzconj")
            local_ndof = 2*order+1;
        else
            local_ndof = BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1);
        nel = ma->GetNE();
        ndof = local_ndof * nel;

        SetDefinedOn(BND, BitArray(ma->GetNRegions(BND)).Clear());

        // evaluators
        switch (D)
        {
            case 2:
                {
                    if(eqtyp=="fowave" || eqtyp=="foqtwave")
                    {
                        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<2,BlockMappedElement<2>>>>();
                        //flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>>();
                        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
                    }
                    else if(eqtyp=="helmholtz" || eqtyp=="helmholtzconj")
                    {
                        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedComplex<2>>>();
                        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradientComplex<2>>>();
                    }
                    else
                    {
                        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<2>>>();
                        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<2>>>();
                        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
                    }
                    break;
                }
            case 3:
                {
                    if(eqtyp=="fowave" || eqtyp=="foqtwave")
                    {
                        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<3,BlockMappedElement<3>>>>();
                        //flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>>();
                        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
                    }
                    else
                    {
                        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMapped<3>>>();
                        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpMappedGradient<3>>>();
                        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
                    }
                    break;
                }
        }
        // basis
        switch (D)
        {
            case 2:
                {
                    if(eqtyp=="laplace")
                        basismat = TLapBasis<2>::Basis(order, basistype);
                    else if(eqtyp=="fowave" || eqtyp=="foqtwave"){
                        basismats.SetSize(D);
                        for(int d=0;d<D;d++) basismats[d]=FOTWaveBasis<1>::Basis(order, d);
                        basis = new FOQTWaveBasis<1>;
                    }
                    else
                    {
                        basismat = TWaveBasis<1>::Basis(order, basistype,eqtyp=="fowave_reduced");
                        basis = new QTWaveBasis<1>;
                    }
                    break;
                }
            case 3:
                {
                    if(eqtyp=="laplace")
                        basismat = TLapBasis<3>::Basis(order, basistype);
                    else if(eqtyp=="fowave" || eqtyp=="foqtwave"){
                        basismats.SetSize(D);
                        for(int d=0;d<D;d++) basismats[d]=FOTWaveBasis<2>::Basis(order, d);
                        basis = new FOQTWaveBasis<2>;
                    }
                    else
                    {
                        basismat = TWaveBasis<2>::Basis(order, basistype,eqtyp=="fowave_reduced");
                        basis = new QTWaveBasis<2>;
                    }
                    break;
                }
        }
    }

    void TrefftzFESpace :: SetWavespeed(shared_ptr<CoefficientFunction> awavespeedcf, shared_ptr<CoefficientFunction> aBBcf, shared_ptr<CoefficientFunction> aGGcf)
    {
        wavespeedcf=awavespeedcf;
        if(aBBcf || eqtyp=="qtwave" || eqtyp=="foqtwave")
        {
            //wavespeedcf = UnaryOpCF(aBBcf/awavespeedcf,GenericSqrt());
            cout << "started auto diff.... ";
            shared_ptr<CoefficientFunction> GGcf = make_shared<ConstantCoefficientFunction>(1)/(awavespeedcf*awavespeedcf);
            shared_ptr<CoefficientFunction> GGcfx = make_shared<ConstantCoefficientFunction>(1)/(awavespeedcf*awavespeedcf);
            if(aGGcf)
            {
                GGcf = aGGcf;
                GGcfx = aGGcf;
            }

            static Timer timerder("QTrefftzDerivatives");
            static Timer timereval("QTrefftzDerEval");
            timerder.Start();
            GGder.SetSize(this->order-(eqtyp=="qtwave"),pow(this->order-(eqtyp=="qtwavw"),D==3));

            for(int ny=0;ny<=(this->order-1-(eqtyp=="qtwave"))*(D==3);ny++)
            {
                for(int nx=0;nx<=this->order-1-(eqtyp=="qtwave");nx++)
                {
                    GGder(nx,ny) = GGcfx;
                    GGcfx = GGcfx->Diff(MakeCoordinateCoefficientFunction(0).get(), make_shared<ConstantCoefficientFunction>(1) );
                }
                GGcf = GGcf->Diff(MakeCoordinateCoefficientFunction(1).get(), make_shared<ConstantCoefficientFunction>(1) );
                GGcfx = GGcf;
            }
            timerder.Stop();


            if(!aBBcf){
                aBBcf = make_shared<ConstantCoefficientFunction>(1);
                cout << "SETTING BB TO 1" << endl;
            }
            static Timer timerbb("QTrefftzBB");
            timerbb.Start();
            shared_ptr<CoefficientFunction> BBcf = aBBcf;
            shared_ptr<CoefficientFunction> BBcfx = aBBcf;
            BBder.SetSize(this->order,(this->order-1)*(D==3)+1);
            for(int ny=0;ny<=(this->order-1)*(D==3);ny++)
            {
                for(int nx=0;nx<=this->order-1;nx++)
                {
                    BBder(nx,ny) = BBcfx;
                    BBcfx = BBcfx->Diff(MakeCoordinateCoefficientFunction(0).get(), make_shared<ConstantCoefficientFunction>(1) );
                }
                BBcf = BBcf->Diff(MakeCoordinateCoefficientFunction(1).get(), make_shared<ConstantCoefficientFunction>(1) );
                BBcfx = BBcf;
            }
            timerbb.Stop();
            cout << "finish" << endl;
        }
    }


    void TrefftzFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
        dnums.SetSize(0);
        if (!DefinedOn (ei) || ei.VB() != VOL) return;
        for(size_t j = ei.Nr()*local_ndof; j<local_ndof*(ei.Nr()+1); j++)
            dnums.Append (j);
    }


    FiniteElement & TrefftzFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
        Ngs_Element ngel = ma->GetElement(ei);
        ELEMENT_TYPE eltype = ngel.GetType();

        if (ei.IsVolume())
        {
            switch (ma->GetElType(ei)) {
                case ET_POINT:
                case ET_SEGM:
                    {
                        throw Exception("illegal dim for space-time element");
                        break;
                    }
                case ET_QUAD:
                case ET_TRIG:
                    {
                        if(eqtyp=="qtwave")
                        {
                            CSR basismat = static_cast<QTWaveBasis<1>*>(basis)->Basis(order, ElCenter<1>(ei), GGder, BBder);
                            return *(new (alloc) ScalarMappedElement<2>(local_ndof,order,basismat,eltype,ElCenter<1>(ei),1.0));
                        }
                        else if(eqtyp=="foqtwave"){
                            Vec<2,CSR> qbasismats;
                            for(int d=0;d<D;d++)
                                qbasismats[d]=
                                    static_cast<FOQTWaveBasis<1>*>(basis)->Basis(order, d, ElCenter<1>(ei), GGder, BBder);
                            return *(new (alloc) BlockMappedElement<2>(local_ndof,order,qbasismats,eltype,ElCenter<1>(ei),1.0));
                        }
                        else if(eqtyp=="fowave"){
                            return *(new (alloc) BlockMappedElement<2>(local_ndof,order,basismats,eltype,ElCenter<1>(ei),Adiam<1>(ei,c),c));
                        }
                        else if(eqtyp=="helmholtz" || eqtyp=="helmholtzconj")
                        {

                            return *(new (alloc) PlaneWaveElement<2>(local_ndof,order,eltype,ElCenter<1>(ei),Adiam<1>(ei,c),c,(eqtyp=="helmholtz"?1:-1)));
                        }
                        else
                            return *(new (alloc) ScalarMappedElement<2>(local_ndof,order,basismat,eltype,ElCenter<1>(ei),Adiam<1>(ei,c),c));
                        break;
                    }
                case ET_HEX:
                case ET_PRISM:
                case ET_PYRAMID:
                case ET_TET:
                    {
                        if(eqtyp=="qtwave")
                        {
                            CSR basismat = static_cast<QTWaveBasis<2>*>(basis)->Basis(order, ElCenter<2>(ei), GGder, BBder);
                            return *(new (alloc) ScalarMappedElement<3>(local_ndof,order,basismat,eltype,ElCenter<2>(ei),1.0));
                        }
                        else if(eqtyp=="foqtwave"){
                            Vec<3,CSR> qbasismats;
                            for(int d=0;d<D;d++)
                                qbasismats[d]=
                                    static_cast<FOQTWaveBasis<2>*>(basis)->Basis(order, d, ElCenter<2>(ei), GGder, BBder);
                            return *(new (alloc) BlockMappedElement<3>(local_ndof,order,qbasismats,eltype,ElCenter<2>(ei),1.0));
                        }
                        else if(eqtyp=="fowave"){
                            return *(new (alloc) BlockMappedElement<3>(local_ndof,order,basismats,eltype,ElCenter<2>(ei),Adiam<2>(ei,c),c));
                        }
                        else
                            return *(new (alloc) ScalarMappedElement<3>(local_ndof,order,basismat,eltype,ElCenter<2>(ei),Adiam<2>(ei,c),c));
                    }
                    break;
            }
        }
        //else
        try
        {
            return SwitchET<ET_POINT,ET_SEGM,ET_TRIG,ET_QUAD>
                (eltype,
                 [&alloc] (auto et) -> FiniteElement&
                 { return * new (alloc) DummyFE<et.ElementType()>; });
        }
        catch (Exception &e)
        {
            throw Exception("illegal element type in Trefftz::GetSurfaceFE");
        }
    }


    DocInfo TrefftzFESpace :: GetDocu ()
    {
        //auto docu = FESpace::GetDocu();
        DocInfo docu;
        docu.short_docu = "Trefftz space for different PDEs. Use kwarg 'eq' to choose the PDE, currently implemented are:\n"
            " - laplace - for Laplace equation\n"
            " - wave - for the second order acoustic wave equation\n"
            " - fowave - for the first order acoustic wave equation, returns TnT (sigv,tauw)\n"
            " - helmholtz - planewaves for the helmholtz equation\n"
            " - helmholtzconj - returns the complex conjungate of the planewaves \n";
        docu.Arg("eq") = "string\n"
            "  Choose type of Trefftz functions.";
        docu.Arg("order") = "int = 1\n"
          "  Order of finite element space";
        docu.Arg("dgjumps") = "bool = True\n"
          "  Enable discontinuous space for DG methods, this flag is always True for trefftzfespace.";
        docu.Arg("complex") = "bool = False\n"
          "  Set if FESpace should be complex";
        //docu.Arg("useshift") = "bool = True\n"
            //"  use shift of basis functins to element center and scale them";
        //docu.Arg("basistype")
        //docu.Arg("wavespeed")
        return docu;
    }

    static RegisterFESpace<TrefftzFESpace> initi_trefftz ("trefftzfespace");




    //////////////////////////// Trefftz basis ////////////////////////////

    template<int D,typename TFunc>
    void TraversePol( int order, const TFunc &func)
    {
        switch(D)
        {
            case 0:
                break;
            case 1:
                for (int i = 0, ii = 0; i <= order; i++)
                                func(ii++,Vec<1,int>{i});
                break;
            case 2:
                for (int i = 0, ii = 0; i <= order; i++)
                    for (int j = 0; j <= order-i; j++)
                                func(ii++,Vec<2,int>{j,i});
                break;
            case 3:
                for (int i = 0, ii = 0; i <= order; i++)
                    for (int j = 0; j <= order-i; j++)
                        for (int k = 0; k <= order-i-j; k++)
                                func(ii++,Vec<3,int>{k,j,i});
                break;
            case 4:
                for (int i = 0, ii = 0; i <= order; i++)
                    for (int j = 0; j <= order-i; j++)
                        for (int k = 0; k <= order-i-j; k++)
                            for (int l = 0; l <= order-i-j-k; l++)
                                func(ii++,Vec<4,int>{l,k,j,i});
                break;
            default:
                throw Exception("TraverseDimensions: too many dimensions!");
        }
    }

    // k-th coeff of Legendre polynomial of degree n in monomial basis
    constexpr double LegCoeffMonBasis(int n, int k)
    {
        if(n==0) return 1;
        if(k>n) return 0;
        if((n+k)%2) return 0;
        double coeff = pow(2,-n) * pow(-1,floor((n-k)/2)) * BinCoeff(n,floor((n-k)/2)) * BinCoeff(n+k,n);
        // double coeff = pow(2,-n) * pow(-1,k) * BinCoeff(n,k) * BinCoeff(2*n-2*k,n);
        return coeff;
    }

    // k-th coeff of Chebyshev polynomial of degree n in monomial basis
    constexpr double ChebCoeffMonBasis(int n, int k)
    {
        if(n==0) return 1;
        if(k>n) return 0;
        if((n+k)%2) return 0;
        double coeff = pow(2,k-1)*n*pow(-1,floor((n-k)/2)) * tgamma((n+k)/2)/(tgamma(floor((n-k)/2)+1)*tgamma(k+1));
        return coeff;
    }

    constexpr int factorial(int n)
    {
        return n>1 ? n * factorial(n-1) : 1;
    }



    template<int D>
    CSR TWaveBasis<D> :: Basis(int ord, int basistype, int fowave)
    {
        CSR tb;
        const int ndof = (BinCoeff(D + ord, ord) + BinCoeff(D + ord-1, ord-1));
        const int npoly = (BinCoeff(D+1 + ord, ord));
        Matrix<> trefftzbasis(ndof,npoly);
        trefftzbasis = 0;
        for(int basis=0;basis<ndof;basis++)
        {
            int tracker = 0;
            TraversePol<D+1>( ord, [&](int i, Vec<D+1,int> coeff){
                if(tracker >= 0) tracker++;
                int indexmap = PolBasis::IndexMap2<D>(coeff, ord);
                int k = coeff(D);
                if(k==0 || k==1)
                {
                    switch (basistype) {
                        case 0:
                            if(tracker>basis)
                            {
                                //trefftzbasis( i, setbasis++ ) = 1.0; //set the l-th coeff to 1
                                trefftzbasis(basis,indexmap) = 1;
                                tracker = -1;
                            }
                            //i += ndof-1;	//jump to time = 2 if i=0
                            break;
                        case 1:
                            if((k == 0 && basis < BinCoeff(D + ord, ord)) || (k == 1 && basis >= BinCoeff(D + ord, ord))){
                                trefftzbasis( basis,indexmap ) = 1;
                                for(int exponent :  coeff.Range(0,D)) trefftzbasis( basis,indexmap ) *= LegCoeffMonBasis(basis,exponent);}
                            break;
                        case 2:
                            if((k == 0 && basis < BinCoeff(D + ord, ord)) || (k == 1 && basis >= BinCoeff(D + ord, ord))){
                                trefftzbasis( basis,indexmap ) = 1;
                                for(int exponent :  coeff.Range(0,D)) trefftzbasis( basis,indexmap ) *= ChebCoeffMonBasis(basis,exponent);}
                            break;
                    }
                }
                else if(coeff(D)>1)
                {
                    for(int m=0;m<D;m++) //rekursive sum
                    {
                        Vec<D+1, int> get_coeff = coeff;
                        get_coeff[D] = get_coeff[D] - 2;
                        get_coeff[m] = get_coeff[m] + 2;
                        trefftzbasis( basis, indexmap) += (coeff(m)+1) * (coeff(m)+2) * trefftzbasis(basis, PolBasis::IndexMap2<D>(get_coeff, ord));
                    }
                    double wavespeed = 1.0;
                    trefftzbasis(basis, indexmap) *= wavespeed*wavespeed/(k * (k-1));
                }
            });
        }
        MatToCSR(trefftzbasis.Rows(fowave,ndof),tb);
        return tb;
    }

    template class TWaveBasis<1>;
    template class TWaveBasis<2>;
    template class TWaveBasis<3>;


    template<int D>
    CSR TLapBasis<D> :: Basis(int ord, int basistype)
    {
        CSR tb;
        const int ndof = (BinCoeff(D-1 + ord, ord) + BinCoeff(D-1 + ord-1, ord-1));
        const int npoly = (BinCoeff(D + ord, ord));
        Matrix<> trefftzbasis(ndof,npoly);
        trefftzbasis = 0;
        for(int basis=0;basis<ndof;basis++)
        {
            int tracker = 0;
            TraversePol<D>( ord, [&](int i, Vec<D,int> coeff){
                if(tracker >= 0) tracker++;
                int indexmap = PolBasis::IndexMap2<D-1>(coeff, ord);
                int k = coeff(D-1);
                if(k==0 || k==1)
                {
                    if(tracker>basis)
                    {
                        //trefftzbasis( i, setbasis++ ) = 1.0; //set the l-th coeff to 1
                        trefftzbasis(basis,indexmap) = 1;
                        tracker = -1;
                    }
                }
                else if(coeff(D-1)>1)
                {
                    for(int m=0;m<D-1;m++) //rekursive sum
                    {
                        Vec<D, int> get_coeff = coeff;
                        get_coeff[D-1] = get_coeff[D-1] - 2;
                        get_coeff[m] = get_coeff[m] + 2;
                        trefftzbasis( basis, indexmap) -= (coeff(m)+1) * (coeff(m)+2) * trefftzbasis(basis, PolBasis::IndexMap2<D-1>(get_coeff, ord));
                    }
                    double lapcoeff = 1.0;
                    trefftzbasis(basis, indexmap) *= lapcoeff*lapcoeff/(k * (k-1));
                }
            });
        }
        MatToCSR(trefftzbasis,tb);
        return tb;
    }

    template class TLapBasis<1>;
    template class TLapBasis<2>;
    template class TLapBasis<3>;


    template<int D>
    CSR FOTWaveBasis<D> :: Basis(int ord, int rdim)
    {
        const int ndof = (D+1) * BinCoeff(ord + D, D);
        const int npoly = BinCoeff((D+1) + ord, ord);
        Array<Matrix<>> trefftzbasis(D+1);
        for (int d = 0; d < D+1; d++)
        {
            trefftzbasis[d].SetSize(ndof,npoly);
            trefftzbasis[d] = 0;
        }
        for(int basis=0;basis<ndof;basis++)
        {
            int tracker = 0;
            TraversePol<D+1>( ord, [&](int i, Vec<D+1,int> coeff){
                if(tracker >= 0) tracker++;
                int indexmap = PolBasis::IndexMap2<D>(coeff, ord);
                if(coeff(D)==0 && tracker*(D+1)>basis)
                {
                    int d = basis % (D+1);
                    trefftzbasis[d](basis,indexmap) = 1;
                    tracker = -1;
                }
                else if(coeff(D)>0)
                {
                    int k = coeff(D);
                    for(int d=0;d<D;d++)
                    {
                        Vec<D+1, int> get_coeff = coeff;
                        get_coeff[D] = get_coeff[D] - 1;
                        get_coeff[d] = get_coeff[d] + 1;
                        trefftzbasis[d](basis, indexmap) = (-1.0/k)*(coeff(d)+1) * trefftzbasis[D](basis, PolBasis::IndexMap2<D>(get_coeff, ord));
                        trefftzbasis[D]( basis, indexmap) += (-1.0/k)*(coeff(d)+1) * trefftzbasis[d](basis, PolBasis::IndexMap2<D>(get_coeff, ord));
                    }
                }

            });
        }
        Array<CSR> tb(D+1);
        for (int d = 0; d < D+1; d++)
        {
            //cout << d << endl << trefftzbasis[d] << endl;
            MatToCSR(trefftzbasis[d],tb[d]);
        }
        return tb[rdim];
    }

    template class FOTWaveBasis<1>;
    template class FOTWaveBasis<2>;
    template class FOTWaveBasis<3>;


    //////////////////////////// quasi-Trefftz basis ////////////////////////////

    template<int D>
    CSR QTWaveBasis<D> :: Basis(int ord, Vec<D+1> ElCenter, Matrix<shared_ptr<CoefficientFunction>> GGder, Matrix<shared_ptr<CoefficientFunction>> BBder, double elsize, int basistype)
    {
        lock_guard<mutex> lock(gentrefftzbasis);
        string encode = to_string(ord) + to_string(elsize);
        for(int i=0;i<D;i++)
            encode += to_string(ElCenter[i]);

        if ( gtbstore[encode][0].Size() == 0)
        {
            IntegrationPoint ip(ElCenter,0);
            Mat<D,D> dummy;
            FE_ElementTransformation<D,D> et(D==3?ET_TET:D==2?ET_TRIG:ET_SEGM,dummy);
            MappedIntegrationPoint<D,D> mip(ip,et,0);
            for(int i=0;i<D;i++)
                mip.Point()[i] = ElCenter[i];

            Matrix<> BB(ord,(ord-1)*(D==2)+1);
            Matrix<> GG(ord-1,(ord-2)*(D==2)+1);
            for(int ny=0;ny<=(ord-1)*(D==2);ny++)
            {
                for(int nx=0;nx<=ord-1;nx++)
                {
                    double fac = (factorial(nx)*factorial(ny));
                    BB(nx,ny) = BBder(nx,ny)->Evaluate(mip)/fac * pow(elsize,nx+ny);
                    if(nx<ord-1 && ny<ord-1)
                        GG(nx,ny) = GGder(nx,ny)->Evaluate(mip)/fac * pow(elsize,nx+ny);
                }
            }

            const int ndof = (BinCoeff(D + ord, ord) + BinCoeff(D + ord-1, ord-1));
            const int npoly = BinCoeff(D+1 + ord, ord);
            Matrix<> qbasis(ndof,npoly);
            qbasis = 0;

            for(int t=0, basisn=0;t<2;t++)
                for(int x=0;x<=ord-t;x++)
                    for(int y=0;y<=(ord-x-t)*(D==2);y++)
                    {
                        Vec<D+1, int> index;
                        index[D] = t;
                        index[0] = x;
                        if(D==2) index[1]=y;
                        qbasis( basisn++, PolBasis::IndexMap2<D>(index, ord))=1.0;
                    }

            for(int basisn=0;basisn<ndof;basisn++)
            {
                for(int ell=0;ell<ord-1;ell++)
                {
                    for(int t=0;t<=ell;t++)
                    {
                        for(int x=(D==1?ell-t:0);x<=ell-t;x++)
                        {
                            int y = ell-t-x;
                            Vec<D+1, int> index;
                            index[1] = y; index[0] = x; index[D] = t+2;
                            double* newcoeff =& qbasis( basisn, PolBasis::IndexMap2<D>(index, ord));
                            *newcoeff = 0;

                            for(int betax=0;betax<=x;betax++)
                                for(int betay=(D==2)?0:y;betay<=y;betay++)
                                {
                                    index[1] = betay; index[0] = betax+1; index[D] = t;
                                    int getcoeffx = PolBasis::IndexMap2<D>(index, ord);
                                    index[1] = betay+1; index[0] = betax; index[D] = t;
                                    int getcoeffy = PolBasis::IndexMap2<D>(index, ord);
                                    index[1] = betay; index[0] = betax+2; index[D] = t;
                                    int getcoeffxx = PolBasis::IndexMap2<D>(index, ord);
                                    index[1] = betay+2; index[0] = betax; index[D] = t;
                                    int getcoeffyy = PolBasis::IndexMap2<D>(index, ord);

                                    *newcoeff +=
                                        (betax+2)*(betax+1)/((t+2)*(t+1)*GG(0)) * BB(x-betax,y-betay)
                                        * qbasis( basisn, getcoeffxx)
                                        + (x-betax+1)*(betax+1)/((t+2)*(t+1)*GG(0)) * BB(x-betax+1,y-betay)
                                        * qbasis( basisn, getcoeffx);
                                    if(D==2)
                                        *newcoeff +=
                                            (betay+2)*(betay+1)/((t+2)*(t+1)*GG(0)) * BB(x-betax,y-betay)
                                            * qbasis( basisn, getcoeffyy)
                                            + (y-betay+1)*(betay+1)/((t+2)*(t+1)*GG(0)) * BB(x-betax,y-betay+1)
                                            * qbasis( basisn, getcoeffy);
                                    if(betax+betay == x+y) continue;
                                    index[1] = betay; index[0] = betax; index[D] = t+2;
                                    int getcoeff = PolBasis::IndexMap2<D>(index, ord);

                                    *newcoeff
                                        -= GG(x-betax,y-betay)*qbasis( basisn, getcoeff) / GG(0);
                                }
                        }
                    }
                }
            }

            MatToCSR(qbasis,gtbstore[encode]);
        }

        if ( gtbstore[encode].Size() == 0)
        {
            stringstream str;
            str << "failed to generate trefftz basis of order " << ord << endl;
            throw Exception (str.str());
        }

        return gtbstore[encode];
    }

    template class QTWaveBasis<1>;
    template class QTWaveBasis<2>;


    template<int D>
    CSR FOQTWaveBasis<D> :: Basis(int ord, int rdim, Vec<D+1> ElCenter, Matrix<shared_ptr<CoefficientFunction>> GGder, Matrix<shared_ptr<CoefficientFunction>> BBder, double elsize)
    {
        lock_guard<mutex> lock(gentrefftzbasis);
        string encode = to_string(ord) + to_string(elsize);
        for(int i=0;i<D;i++)
            encode += to_string(ElCenter[i]);

        if ( gtbstore[0][encode][0].Size() == 0)
        {
            IntegrationPoint ip(ElCenter,0);
            Mat<D,D> dummy;
            FE_ElementTransformation<D,D> et(D==3?ET_TET:D==2?ET_TRIG:ET_SEGM,dummy);
            MappedIntegrationPoint<D,D> mip(ip,et,0);
            for(int i=0;i<D;i++)
                mip.Point()[i] = ElCenter[i];

            Matrix<> BB(ord,(ord-1)*(D==2)+1);
            Matrix<> GG(ord,(ord-1)*(D==2)+1);
            for(int ny=0;ny<=(ord-1)*(D==2);ny++)
            {
                for(int nx=0;nx<=ord-1;nx++)
                {
                    double fac = (factorial(nx)*factorial(ny));
                    BB(nx,ny) = BBder(nx,ny)->Evaluate(mip)/fac * pow(elsize,nx+ny);
                    GG(nx,ny) = GGder(nx,ny)->Evaluate(mip)/fac * pow(elsize,nx+ny);
                }
            }

            const int ndof = (D+1) * BinCoeff(ord + D, D);
            const int npoly = BinCoeff((D+1) + ord, ord);
            Array<Matrix<>> qbasis(D+1);
            for(int d=0; d<D+1;d++)
            {
                qbasis[d].SetSize(ndof,npoly);
                qbasis[d] = 0;
            }

            for(int d=0, basisn=0; d<D+1;d++)
            {
                for(int x=0;x<=ord;x++)
                    for(int y=0;y<=(ord-x)*(D==2);y++)
                    {
                        Vec<D+1, int> index;
                        index[1]=y;
                        index[0] = x;
                        index[D] = 0;
                        qbasis[d]( basisn++, PolBasis::IndexMap2<D>(index, ord))=1.0;
                    }
            }


            for(int basisn=0;basisn<ndof;basisn++)
            {
                for(int ell=0;ell<ord;ell++)
                {
                    for(int t=0;t<=ell;t++)
                    {
                        for(int x=(D==1?ell-t:0);x<=ell-t;x++)
                        {
                            int y = ell-t-x;
                            Vec<D+1, int> index;
                            index[1] = y; index[0] = x; index[D] = t+1;
                            int newindex = PolBasis::IndexMap2<D>(index, ord);
                            double* newcoefft =& qbasis[D]( basisn, newindex);
                            for(int d=0; d<D;d++)
                            {
                                double* newcoeff =& qbasis[d]( basisn, newindex);

                                index[1] = y+(d==1); index[0] = x+(d==0); index[D] = t;
                                int getcoeff = PolBasis::IndexMap2<D>(index, ord);
                                *newcoeff
                                    = -qbasis[D]( basisn, getcoeff)*index[d]/(t+1)/BB(0);
                                *newcoefft
                                    -= qbasis[d]( basisn, getcoeff)*index[d]/(t+1)/GG(0);
                                for(int betax=0;betax<=x;betax++)
                                    for(int betay=(D==1)?y:0;betay<=y;betay++)
                                    {
                                        if(betax+betay==x+y) continue;
                                        index[1] = betay; index[0] = betax; index[D] = t+1;
                                        int getcoeff = PolBasis::IndexMap2<D>(index, ord);
                                        *newcoeff
                                            -= BB(x-betax,y-betay)*qbasis[d]( basisn, getcoeff) / BB(0);
                                        if(d==0)
                                        *newcoefft
                                            -= GG(x-betax,y-betay)*qbasis[D]( basisn, getcoeff) / GG(0);
                                    }
                            }
                        }
                    }
                }
            }

            for (int d = 0; d < D+1; d++)
            {
                MatToCSR(qbasis[d],gtbstore[d][encode]);
            }
        }

        if ( gtbstore[0][encode].Size() == 0)
        {
            stringstream str;
            str << "failed to generate trefftz basis of order " << ord << endl;
            throw Exception (str.str());
        }

        return gtbstore[rdim][encode];

    }

    template class FOQTWaveBasis<1>;
    template class FOQTWaveBasis<2>;

}


#ifdef NGS_PYTHON
void ExportTrefftzFESpace(py::module m)
{
    using namespace ngcomp;
    //using namespace ngfem;
    //[>
    //We just export the class here and use the FESpace constructor to create our space.
    //This has the advantage, that we do not need to specify all the flags to parse (like
    //dirichlet, definedon,...), but we can still append new functions only for that space.
    //*/
    //py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>, FESpace>
    //(m, "TrefftzFESpace", "FESpace with first order and second order trigs on 2d mesh")
    //.def("GetNDof", &TrefftzFESpace::GetNDof)
    //;
    //m.def("GetNDof", [](shared_ptr<FESpace> fes) {
    //cout << typeid(*fes).name() << endl;
    ////fes->GetNDof();
    //});

    ExportFESpace<TrefftzFESpace>(m, "trefftzfespace")
        .def("GetDocu", &TrefftzFESpace::GetDocu)
        .def("GetNDof", &TrefftzFESpace::GetNDof)
        .def("SetWavespeed", &TrefftzFESpace::SetWavespeed, py::arg("Wavespeed"), py::arg("BBcf")=nullptr, py::arg("GGcf")=nullptr)
        ;

   //ExportFESpace<FOTWaveFESpace, CompoundFESpace> (m, "FOTWave");
}
#endif // NGS_PYTHON
