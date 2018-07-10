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
    Vec<D+1> VertexTimesTentFace(Tent tent, const INT<D+1>& verts )
    {
        // determine linear basis function to use for tent face
        Vec<D+1> bs=0;
        for(int ivert = 0;ivert<verts.Size();ivert++)
        {
            if(verts[ivert]==tent.vertex) bs[ivert]=tent.ttop;
            for (int k = 0; k < tent.nbv.Size(); k++)
                if(verts[ivert] == tent.nbv[k]) bs[ivert] =  tent.nbtime[k];
        }
        return bs;
    }

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
            cout << "vertex: " << tent.vertex << ", tbot = " << tent.tbot << ", ttop = " << tent.ttop << endl;
            cout << "neighbour vertices: " << endl;
            for (int k = 0; k < tent.nbv.Size(); k++)
                cout << k << ": " << tent.nbv[k] << " " << tent.nbtime[k] << endl;
            //cout << "elements: " << endl << tent.els << endl;

            int nbasis = BinCoeff(D-1 + order, order) + BinCoeff(D-1 + order-1, order-1);
            Matrix<double> elmatrix(nbasis,nbasis);
            Vector<double> elvector(nbasis);

            for(auto el : tent.els)
            {
                MappedIntegrationRule<1,D> mir(ir, ma->GetTrafo(el,lh), lh); // <dim  el, dim space>

                INT<D+1> verts = ma->GetEdgePNums(el);
                Vec<D+1> bs = VertexTimesTentFace<D>(tent, verts);

                auto p = ma->GetPoint<1>(tent.vertex);
                cout << "ns: " << bs << endl;
                //cout << "mir0: " << mir[0] << " val: " << faceint.Evaluate( ir[0], bs) << endl;
                //cout << "mir1: " << mir[1] << " val: " << faceint.Evaluate( ir[1], bs) << endl;
                //cout << "mir1: " << mir[ir.Size()-1]<<mir[ir.Size()-1].GetWeight() << endl;

                for(int imip=0;imip<mir.Size();imip++)
                {
                    auto p = mir[imip].GetPoint();
                    Vec<D+1> tmip;                     
                    tmip.Range(0,D) = p; 
                    tmip(D) = faceint.Evaluate(ir[imip], bs);
                    //tent.nbtime[]
                    cout << "tmip: " << tmip << endl;
                }
                cout << "edge verts: " << ma->GetPoint<1>(verts[0]) << ma->GetPoint<1>(verts[1]) << endl; 
            }

        });
        // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
        // py::object ffes = py::cast(fes);
        // auto pyspace = py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str());
        // py::object pyfes = et.attr("GetFESTrefftz")(ma);
        // FESpace *ffes = pyfes.cast<FESpace *>();
        // et.attr("EvolveTent")(pyfes,?,?);
    }
    /*
       template<int D>
       double AreaTentFace(Tent tent, int elnr, shared_ptr<MeshAccess> ma) 
       {
       Mat<D+2,D+2> volmat;
       volmat.Row(0) = 1;
       volmat.Col(0) = 1;
       volmat(0,0) = 0;
       auto verts = ma->GetEdgePNums(elnr);
       Vec<D+1> bs =  VertexTimesTentFace<D>(tent, verts);

       for(int i=1;i<verts.Size();i++)
       {

       }
       double bla = 3.0;
       return bla;

       double anisotropicdiam = 0.0;
       auto vertices_index = ma->GetElVertices(ei);
       for(auto vertex1 : vertices_index)
       {
       for(auto vertex2 : vertices_index)
       {
       Vec<D> v1 = ma->GetPoint<D>(vertex1);
       Vec<D> v2 = ma->GetPoint<D>(vertex2);
    //cout << "v1: " << v1 << " v1 part: " << v1(1,D-1) << "norm " << L2Norm(v1) << endl ;
    anisotropicdiam = max( anisotropicdiam, sqrt( L2Norm2(v1(0,D-2) - v2(0,D-2)) + pow(c*(v1(D-1)-v2(D-1)),2) ) );
    }
    }
    return anisotropicdiam;
    }
    */
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
              Vector<double> wavefront(ir_order * ne * (D+2));
              //return EvolveTents<1>(ma,wavespeed,dt,wavefront);
          }//, py::call_guard<py::gil_scoped_release>()
         );
}
#endif // NGS_PYTHON

