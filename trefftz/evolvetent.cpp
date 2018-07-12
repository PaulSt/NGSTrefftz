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


    template<int D>    
    void EvolveTents(shared_ptr<MeshAccess> ma, double wavespeed, double dt, Vector<double> wavefront)
    {
        int order = 3;
        int ir_order = ceil((order+1)/1);
        int ne = ma->GetNE();

        // wavefront.Resize(ir_order * ma->GetNEdges() * (D+2));
        T_TrefftzElement<D+1> tel(order,wavespeed);
        IntegrationRule ir(ET_SEGM, ir_order);
        ScalarFE<ET_SEGM,D> faceint;

        //py::module et = py::module::import("DGeq");

        TentPitchedSlab<D> tps = TentPitchedSlab<D>(ma);      // collection of tents in timeslab
        tps.PitchTents(dt, wavespeed); // adt = time slab height, wavespeed
        LocalHeap lh(order * D * 1000);

        RunParallelDependency (tps.tent_dependency, [&] (int i) {
            HeapReset hr(lh);
            // LocalHeap slh = lh.Split();  // split to threads
            Tent & tent = *tps.tents[i];
            //cout << tent << endl;
            cout << endl << "tent: " << i << " " << tent.vertex << endl;
            cout << "vertex: " << tent.vertex << " at: " << ma->GetPoint<D>(tent.vertex) <<  ", tbot = " << tent.tbot << ", ttop = " << tent.ttop << endl;
            cout << "neighbour vertices: " << endl;
            for (int k = 0; k < tent.nbv.Size(); k++)
                cout << k << ": " << tent.nbv[k] << " at: " << ma->GetPoint<D>(tent.nbv[k]) <<" t: " << tent.nbtime[k] << endl;

            int nbasis = BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1);
            Matrix<double> elmatrix(nbasis,nbasis);
            Vector<double> elvector(nbasis);

            for(auto elnr: tent.els)
            {
                MappedIntegrationRule<1,D> mir(ir, ma->GetTrafo(elnr,lh), lh); // <dim  el, dim space>

                INT<D+1> verts = ma->GetEdgePNums(elnr);
                Vec<D+1> bs = VertexTimesTentFace<D>(tent, verts);

                AreaTentFace<D>(tent, elnr, ma);

                for(int imip=0;imip<mir.Size();imip++)
                {
                    Vec<D+1> tmip_point;
                    tmip_point.Range(0,D) = mir[imip].GetPoint(); 
                    tmip_point(D) = faceint.Evaluate(ir[imip], bs);
                    IntegrationPoint ip(tmip_point, ir[imip].Weight());
                    MappedIntegrationPoint<D+1,D+1> tmip(ip, FE_ElementTransformation<D+1,D+1>() );

                    Matrix<> blaa(D+1,nbasis);
                    tel.CalcDShape(tmip,blaa);

                    cout << "tmip_point: " << tmip_point << endl;
                    
                }
            }

        });
        // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
        // py::object ffes = py::cast(fes);
        // auto pyspace = py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str());
        // py::object pyfes = et.attr("GetFESTrefftz")(ma);
        // FESpace *ffes = pyfes.cast<FESpace *>();
        // et.attr("EvolveTent")(pyfes,?,?);
    }

    template<int D>
    Vec<D+1> VertexTimesTentFace(const Tent& tent, const INT<D+1>& verts)
    {
        Vec<D+1> bs;
        // determine linear basis function coeffs to use for tent face
        for(int ivert = 0;ivert<verts.Size();ivert++)
        {
            if(verts[ivert] == tent.vertex) bs[ivert] = tent.ttop;
            for (int k = 0; k < tent.nbv.Size(); k++)
                if(verts[ivert] == tent.nbv[k]) bs[ivert] = tent.nbtime[k];
        }
        return bs;
    }

    template<int D>
    double AreaTentFace(const Tent& tent, int elnr, shared_ptr<MeshAccess> ma) 
    {
        INT<D+1> verts = ma->GetEdgePNums(elnr);
        Vec<D+1> bs =  VertexTimesTentFace<D>(tent, verts);
        Vec<D+1, Vec<D+1>> v;
        for(int i=0;i<=D;i++)
        {
            v[i].Range(0,D-1) = ma->GetPoint<D>(verts[i]);
            v[i](D) = bs(i);
        }
        switch(D) {
            case 1 : return L2Norm(v[0]-v[1]);
                     break;    
            case 2 :{
                        double a = L2Norm2(v[0]-v[1]);
                        double b = L2Norm2(v[1]-v[2]);
                        double c = L2Norm2(v[0]-v[2]);
                        double s = 0.5*(a+b+c);
                        return sqrt(s*(s-a)*(s-b)*(s-c));
                        break;
                    }
            default: return 0;
        }
    }
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m)
{
    m.def("EvolveTents", [](shared_ptr<MeshAccess> ma, double wavespeed, double dt) //-> shared_ptr<MeshAccess>
          {
              int D=1;
              int order = 3;
              int ir_order = ceil((order+1)/1);
              int ne = ma->GetNE();
              Vector<double> wavefront; //(ir_order * ne * (D+2));
              EvolveTents<1>(ma,wavespeed,dt,wavefront);
          }//, py::call_guard<py::gil_scoped_release>()
         );
}
#endif // NGS_PYTHON

