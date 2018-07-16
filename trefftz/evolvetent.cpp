#include "evolvetent.hpp"
#include "trefftzelement.hpp"
#include "tents/tents.hpp"
#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>

namespace ngcomp
{
  template <int D>
  void EvolveTents (shared_ptr<MeshAccess> ma, double wavespeed, double dt,
                    Vector<double> wavefront)
  {
    int order = 3;
    int ir_order = ceil ((order + 1) / 1);
    int ne = ma->GetNE ();

    // wavefront.Resize(ir_order * ma->GetNEdges() * (D+2));
    T_TrefftzElement<D + 1> tel (order, wavespeed);
    int nbasis = tel.GetNBasis ();
    cout << "NBASIS: " << nbasis << endl;

    IntegrationRule ir (ET_SEGM, ir_order);
    ScalarFE<ET_SEGM, D> faceint;

    // py::module et = py::module::import("DGeq");

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (ma);  // collection of tents in timeslab
    tps.PitchTents (dt, wavespeed); // adt = time slab height, wavespeed
    LocalHeap lh (order * D * 1000);

    RunParallelDependency (tps.tent_dependency, [&] (int i) {
      HeapReset hr (lh);
      // LocalHeap slh = lh.Split();  // split to threads
      Tent *tent = tps.tents[i];
      // cout << tent << endl;
      // cout << endl << "tent: " << i << " " << tent.vertex << endl;
      // cout << "vertex: " << tent.vertex << " at: " <<
      // ma->GetPoint<D>(tent.vertex) <<  ", tbot = " << tent.tbot << ", ttop =
      // " << tent.ttop << endl; cout << "neighbour vertices: " << endl; for
      // (int k = 0; k < tent.nbv.Size(); k++) cout << k << ": " << tent.nbv[k]
      // << " at: " << ma->GetPoint<D>(tent.nbv[k]) <<" t: " << tent.nbtime[k]
      // << endl;

      for (auto elnr : tent->els)
        {
          Matrix<double> elmatrix (nbasis, nbasis);
          Vector<double> elvector (nbasis);
          elmatrix = 0;
          elvector = 0;

          MappedIntegrationRule<1, D> mir (ir, ma->GetTrafo (elnr, lh),
                                           lh); // <dim  el, dim space>

          INT<D + 1> verts = ma->GetEdgePNums (elnr);
          Vec<D + 1> bs = TentFaceVertexTimes<D> (tent, verts);

          Mat<D + 1, D + 1> v = TentFaceVerts<D> (tent, elnr, ma);
          double A = TentFaceArea<D> (v);
          Vec<D + 1> n = TentFaceNormal<D> (v, 1);

          for (int imip = 0; imip < mir.Size (); imip++)
            {
              Vec<D + 1> p;
              p.Range (0, D) = mir[imip].GetPoint ();
              p (D) = faceint.Evaluate (ir[imip], bs);

              Matrix<> dshape (nbasis, D + 1);
              tel.CalcDShape (p, dshape);

              // cout << "A " << A << endl;
              // cout << "n" << n << endl << endl;

              for (int i = 0; i < nbasis; i++)
                {
                  for (int j = 0; j < nbasis; j++)
                    {
                      elmatrix (i, j)
                          += (dshape (i, D) * dshape (j, D) * n (D))
                             * (1 / (wavespeed * wavespeed)) * A
                             * ir[imip].Weight ();
                      elmatrix (i, j)
                          += (InnerProduct (dshape.Row (i).Range (0, D),
                                            dshape.Row (j).Range (0, D))
                              * n (D))
                             * A * ir[imip].Weight ();
                      elmatrix (i, j)
                          += (dshape (i, D)
                              * InnerProduct (dshape.Row (j).Range (0, D),
                                              n.Range (0, D)))
                             * A * ir[imip].Weight ();
                      elmatrix (i, j)
                          += (dshape (j, D)
                              * InnerProduct (dshape.Row (i).Range (0, D),
                                              n.Range (0, D)))
                             * A * ir[imip].Weight ();
                    }
                }

              int offset = elnr * ir_order * (D + 2) + imip * (D + 2);
              for (int j = 0; j < nbasis; j++)
                {
                  elvector (j)
                      -= (wavefront (offset + D + 1) * dshape (j, D) * n (D))
                         * (1 / (wavespeed * wavespeed)) * A
                         * ir[imip].Weight ();
                  elvector (j)
                      -= (InnerProduct (
                              wavefront.Range (offset + 1, offset + D + 1),
                              dshape.Row (j).Range (0, D))
                          * n (D))
                         * A * ir[imip].Weight ();
                  elvector (j) -= (wavefront (offset + D + 1)
                                   * InnerProduct (dshape.Row (j).Range (0, D),
                                                   n.Range (0, D)))
                                  * A * ir[imip].Weight ();
                  elvector (j)
                      -= (dshape (j, D)
                          * InnerProduct (
                              wavefront.Range (offset + 1, offset + D + 1),
                              n.Range (0, D)))
                         * A * ir[imip].Weight ();
                }
            }
          // cout << elmatrix << endl << elvector << endl;
        }
    });
    // std::shared_ptr<FESpace> p = std::make_shared<TrefftzFESpace>(fes);
    // py::object ffes = py::cast(fes);
    // auto pyspace = py::class_<TrefftzFESpace,
    // shared_ptr<TrefftzFESpace>,FESpace> (m, pyname.c_str()); py::object
    // pyfes = et.attr("GetFESTrefftz")(ma); FESpace *ffes = pyfes.cast<FESpace
    // *>(); et.attr("EvolveTent")(pyfes,?,?);
  }

  template <int D>
  Vec<D + 1> TentFaceVertexTimes (Tent *tent, const INT<D + 1> &verts)
  {
    Vec<D + 1> bs;
    // determine linear basis function coeffs to use for tent face
    for (int ivert = 0; ivert < verts.Size (); ivert++)
      {
        if (verts[ivert] == tent->vertex)
          bs[ivert] = tent->ttop;
        for (int k = 0; k < tent->nbv.Size (); k++)
          if (verts[ivert] == tent->nbv[k])
            bs[ivert] = tent->nbtime[k];
      }
    return bs;
  }

  template <int D>
  Mat<D + 1, D + 1>
  TentFaceVerts (Tent *tent, int elnr, shared_ptr<MeshAccess> ma)
  {
    INT<D + 1> verts = ma->GetEdgePNums (elnr);
    Vec<D + 1> bs = TentFaceVertexTimes<D> (tent, verts);
    Mat<D + 1, D + 1> v;
    for (int i = 0; i <= D; i++)
      {
        v.Row (i).Range (0, D) = ma->GetPoint<D> (verts[i]);
        v.Row (i) (D) = bs (i);
      }
    return v;
  }

  template <int D> double TentFaceArea (Mat<D + 1, D + 1> v)
  {
    switch (D)
      {
      case 1:
        return L2Norm (v.Row (0) - v.Row (1));
        break;
      case 2:
        {
          double a = L2Norm2 (v.Row (0) - v.Row (1));
          double b = L2Norm2 (v.Row (1) - v.Row (2));
          double c = L2Norm2 (v.Row (0) - v.Row (2));
          double s = 0.5 * (a + b + c);
          return sqrt (s * (s - a) * (s - b) * (s - c));
          break;
        }
      }
  }

  template <int D> Vec<D + 1> TentFaceNormal (Mat<D + 1, D + 1> v, bool dir)
  {
    Vec<D + 1> normv;
    switch (D)
      {
      case 1:
        {
          normv (0) = v (0, 1) - v (1, 1);
          normv (1) = v (1, 0) - v (0, 0);
          normv /= L2Norm (normv);
          break;
        }
      case 2:
        {
          Vec<D + 1> a = v.Row (0) - v.Row (1);
          Vec<D + 1> b = v.Row (0) - v.Row (2);
          normv (0) = a (1) * b (2) - a (2) * b (1);
          normv (1) = a (2) * b (0) - a (0) * b (2);
          normv (2) = a (0) * b (1) - a (1) * b (0);
          normv /= sqrt (L2Norm2 (a) * L2Norm2 (b)
                         - (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
                               * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]));
          break;
        }
      }
    if (dir == 1)
      normv *= sgn_nozero<double> (normv[D]);
    else
      normv *= (-sgn_nozero<double> (normv[D]));
    return normv;
  }
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m)
{
  m.def ("EvolveTents",
         [] (shared_ptr<MeshAccess> ma, double wavespeed,
             double dt) //-> shared_ptr<MeshAccess>
         {
           int D = 1;
           int order = 3;
           int ir_order = ceil ((order + 1) / 1);
           int ne = ma->GetNE ();
           Vector<double> wavefront (ir_order * ne * (D + 2));
           EvolveTents<1> (ma, wavespeed, dt, wavefront);
         } //, py::call_guard<py::gil_scoped_release>()
  );
}
#endif // NGS_PYTHON
