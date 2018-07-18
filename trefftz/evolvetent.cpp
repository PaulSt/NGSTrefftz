#include "evolvetent.hpp"
#include "trefftzelement.hpp"
#include "tents/tents.hpp"
#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>


namespace ngcomp
{
    template<int D, typename TFUNC>
    Vector<> MakeIC(IntegrationRule ir, shared_ptr<MeshAccess> ma, LocalHeap& lh, TFUNC func){
        Vector<> ic(ir.Size() * ma->GetNE() * (D+2));
        for(int elnr=0;elnr<ma->GetNE();elnr++){
            MappedIntegrationRule<1,D> mir(ir, ma->GetTrafo(elnr,lh), lh); // <dim  el, dim space>
            for(int imip=0;imip<mir.Size();imip++)
            {
                Vec<D> p = mir[imip].GetPoint();
                int offset = elnr*ir.Size()*(D+2) + imip*(D+2);
                ic.Range(offset,offset+D+2) = func(p);
            }
        }
        return ic;
    }

    template<int D>
    void EvolveTents(int order, shared_ptr<MeshAccess> ma, double wavespeed, double dt)
    {
        int ne = ma->GetNE();

        T_TrefftzElement<D+1> tel(order,wavespeed);
        int nbasis = tel.GetNBasis();
        cout << "NBASIS: " << nbasis << endl;

        IntegrationRule ir(ET_SEGM, order);
        int ir_order = ir.Size();//ceil((order+1)/1);
        cout << "irsize: " << ir.Size() << endl;
        ScalarFE<ET_SEGM,D> faceint;

        py::module et = py::module::import("DGeq");

        TentPitchedSlab<D> tps = TentPitchedSlab<D>(ma);      // collection of tents in timeslab
        tps.PitchTents(dt, wavespeed); // adt = time slab height, wavespeed
        LocalHeap lh(order * D * 10000);

        cout << "NE " << ne << " nedg " << ma->GetNEdges() << endl;
        Vector<> wavefront(ir_order * ne * (D+2));
        wavefront = MakeIC<D>(ir,ma,lh,[&](Vec<D> p){
            double x = p[0]; double y = 0;
            Vec<D+2> sol;
            int k = 30;
            sol[1] = (-2*k*((x-0.5)-wavespeed*y));
            sol[2] = (wavespeed*2*k*((x-0.5)-wavespeed*y));
            sol *= exp(-k*((x-0.5)-wavespeed*y)*((x-0.5)-wavespeed*y));
            sol[0] = exp(-k*((x-0.5)-wavespeed*y)*((x-0.5)-wavespeed*y));
            return sol;
        });
        cout << wavefront << endl;

        RunParallelDependency (tps.tent_dependency, [&] (int i) {
            HeapReset hr(lh);
            // LocalHeap slh = lh.Split();  // split to threads
            Tent* tent = tps.tents[i];
            cout << endl << "tent: " << i << " vert: " << tent->vertex << " els: " << tent->els << endl;
            //cout << tent << endl;

            Matrix<double> elmat(nbasis,nbasis);
            Vector<double> elvec(nbasis);
            elmat = 0;
            elvec = 0;
            for(auto elnr: tent->els)
            {
                INT<D+1> vnr = ma->GetEdgePNums(elnr);
                MappedIntegrationRule<1,D> mir(ir, ma->GetTrafo(elnr,lh), lh); // <dim  el, dim space>

                // Integration over top of tent
                cout << "top" << endl;
                Mat<D+1,D+1> v = TentFaceVerts<D>(tent, elnr, ma, 1);
                Vec<D+1> n = TentFaceNormal<D>(v,1);
                Vec<D+1> bs = v.Col(D); 
                double A = TentFaceArea<D>(v);
                for(int imip=0;imip<mir.Size();imip++)
                {
                    Vec<D+1> p;
                    p.Range(0,D) = mir[imip].GetPoint();
                    p(D) = faceint.Evaluate(ir[imip], bs);

                    cout << p << endl;
                    Matrix<> dshape(nbasis,D+1);
                    tel.CalcDShape(p,dshape);

                    for(int i=0;i<nbasis;i++)
                    {
                        for(int j=0;j<nbasis;j++)
                        {
                            elmat(i,j) += ( dshape(i,D)*dshape(j,D)*n(D) ) * (1/(wavespeed*wavespeed)) *A*ir[imip].Weight();
                            elmat(i,j) += ( InnerProduct(dshape.Row(i).Range(0,D),dshape.Row(j).Range(0,D))*n(D) ) *A*ir[imip].Weight();
                            elmat(i,j) += ( dshape(i,D)*InnerProduct(dshape.Row(j).Range(0,D),n.Range(0,D)) ) *A*ir[imip].Weight();
                            elmat(i,j) += ( dshape(j,D)*InnerProduct(dshape.Row(i).Range(0,D),n.Range(0,D)) ) *A*ir[imip].Weight();
                        }
                    }
                }

                // Integration over bot of tent
                cout << "bot" << endl;
                v = TentFaceVerts<D>(tent, elnr, ma, 0);
                n = TentFaceNormal<D>(v,0);
                bs = v.Col(D); 
                A = TentFaceArea<D>(v);
                for(int imip=0;imip<mir.Size();imip++)
                {
                    Vec<D+1> p;
                    p.Range(0,D) = mir[imip].GetPoint();
                    p(D) = faceint.Evaluate(ir[imip], bs);
                    cout << p << endl;

                    Matrix<> dshape(nbasis,D+1);
                    tel.CalcDShape(p,dshape);
                    Vector<> shape(nbasis);
                    tel.CalcShape(p,shape);

                    int offset = elnr*ir_order*(D+2) + imip*(D+2);
                    for(int j=0;j<nbasis;j++)
                    {
                        elvec(j) -= ( wavefront(offset+D+1)*dshape(j,D)*n(D) ) * (1/(wavespeed*wavespeed)) *A*ir[imip].Weight();
                        elvec(j) -= ( InnerProduct(wavefront.Range(offset+1,offset+D+1),dshape.Row(j).Range(0,D))*n(D) ) *A*ir[imip].Weight();
                        elvec(j) -= ( wavefront(offset+D+1)*InnerProduct(dshape.Row(j).Range(0,D),n.Range(0,D)) ) *A*ir[imip].Weight();
                        elvec(j) -= ( dshape(j,D)*InnerProduct(wavefront.Range(offset+1,offset+D+1),n.Range(0,D)) ) *A*ir[imip].Weight();
                        elmat(j) += ( wavefront(offset)*shape(j) ) *A*ir[imip].Weight();

                        for(int i=0;i<nbasis;i++)
                        {
                            elmat(i,j) += ( shape(i)*shape(j) ) *A*ir[imip].Weight();
                        }
                    }
                }
            }

            //cout << elmat << endl << elvec << endl;

            //KrylovSpaceSolver * solver;
            //solver = new CGSolver<double> (elmat);//, pre);
            //solver->SetPrecision(1e-8);
            //solver->SetMaxSteps(200);
            //solver->SetPrintRates (true);
            Matrix<> mat(elmat);
            CalcInverse(elmat);
            Vector<> sol = elmat*elvec;

            for(auto elnr: tent->els)
            {
                INT<D+1> vnr = ma->GetEdgePNums(elnr);
                MappedIntegrationRule<1,D> mir(ir, ma->GetTrafo(elnr,lh), lh); // <dim  el, dim space>

                // eval solution on top of tent
                Mat<D+1,D+1> v = TentFaceVerts<D>(tent, elnr, ma, 1);
                Vec<D+1> n = TentFaceNormal<D>(v,1);
                Vec<D+1> bs = v.Col(D); 
                double A = TentFaceArea<D>(v);
                for(int imip=0;imip<mir.Size();imip++)
                {
                    Vec<D+1> p;
                    p.Range(0,D) = mir[imip].GetPoint();
                    p(D) = faceint.Evaluate(ir[imip], bs);

                    Matrix<> dshape(nbasis,D+1);
                    tel.CalcDShape(p,dshape);
                    Vector<> shape(nbasis);
                    tel.CalcShape(p,shape);

                    int offset = elnr*ir.Size()*(D+2) + imip*(D+2);
                    cout << wavefront(offset) << endl;
                    wavefront(offset) = InnerProduct(shape,sol);
                    cout << wavefront(offset) << endl;
                    wavefront.Range(offset+1,offset+D+2) = dshape*sol;

                }
            }

        });
        cout << wavefront << endl;
        // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
        // py::object ffes = py::cast(fes);
        // auto pyspace = py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str());
        // py::object pyfes = et.attr("GetFESTrefftz")(ma);
        // FESpace *ffes = pyfes.cast<FESpace *>();
        // et.attr("EvolveTent")(pyfes,?,?);
    }

    template<int D>
    Mat<D+1,D+1> TentFaceVerts(Tent* tent, int elnr, shared_ptr<MeshAccess> ma, bool top)
    {
        INT<D+1> vnr = ma->GetEdgePNums(elnr);
        Mat<D+1, D+1> v;
        // determine linear basis function coeffs to use for tent face
        for(int ivert = 0;ivert<vnr.Size();ivert++)
        {
            if(vnr[ivert] == tent->vertex) v(ivert,D) =  top ? tent->ttop : tent->tbot;
            for (int k = 0; k < tent->nbv.Size(); k++)
                if(vnr[ivert] == tent->nbv[k]) v(ivert,D) = tent->nbtime[k];
            v.Row(ivert).Range(0,D) = ma->GetPoint<D>(vnr[ivert]);
        }

        return v;
    }

    template<int D>
    double TentFaceArea( Mat<D+1,D+1> v )
    {
        switch(D) {
            case 1 : return L2Norm(v.Row(0)-v.Row(1));
                     break;
            case 2 :{
                        double a = L2Norm2(v.Row(0)-v.Row(1));
                        double b = L2Norm2(v.Row(1)-v.Row(2));
                        double c = L2Norm2(v.Row(0)-v.Row(2));
                        double s = 0.5*(a+b+c);
                        return sqrt(s*(s-a)*(s-b)*(s-c));
                        break;
                    }
        }
    }

    template<int D>
    Vec<D+1> TentFaceNormal( Mat<D+1,D+1> v, bool top )
    {
        Vec<D+1> normv;
        switch(D){
            case 1: {
                        normv(0) = v(0,1)-v(1,1);
                        normv(1) = v(1,0)-v(0,0);
                        normv /= L2Norm(normv);
                        break;
                    }
            case 2: {
                        Vec<D+1> a = v.Row(0)-v.Row(1);
                        Vec<D+1> b = v.Row(0)-v.Row(2);
                        normv(0) = a(1) * b(2) - a(2) * b(1);
                        normv(1) = a(2) * b(0) - a(0) * b(2);
                        normv(2) = a(0) * b(1) - a(1) * b(0);
                        normv /= sqrt(L2Norm2(a)*L2Norm2(b)- (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]));
                        break;
                    }
        }
        if(top == 1) normv *= sgn_nozero<double>(normv[D]);
        else normv *= (-sgn_nozero<double>(normv[D]));
        return normv;
    }


}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m)
{
    m.def("EvolveTents", [](int order, shared_ptr<MeshAccess> ma, double wavespeed, double dt) //-> shared_ptr<MeshAccess>
          {
              int D=1;
              EvolveTents<1>(order,ma,wavespeed,dt);
          }//, py::call_guard<py::gil_scoped_release>()
         );
}
#endif // NGS_PYTHON

