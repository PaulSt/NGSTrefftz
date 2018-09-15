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

  inline void LapackSolve (SliceMatrix<double> a, SliceVector<double> b)
  {
    integer n = a.Width ();
    integer lda = a.Dist ();
    integer success;
    char trans = 'T';
    integer nrhs = 1;
    ArrayMem<integer, 100> ipiv (n);

    dgetrf_ (&n, &n, &a (0, 0), &lda, &ipiv[0], &success);
    dgetrs_ (&trans, &n, &nrhs, &a (0, 0), &lda, &ipiv[0], &b[0], &lda,
             &success);
    if (success != 0)
      cout << "Lapack error: " << success << endl;
  }

  template <int D>
  void EvolveTents (int order, shared_ptr<MeshAccess> ma, double wavespeed,
                    double dt, SliceVector<double> wavefront, double timeshift)
  {
    LocalHeap lh (100000000);
    T_TrefftzElement<D + 1> tel (order, wavespeed);
    int nbasis = tel.GetNBasis ();

    Vector<> solutioncoeffs (nbasis * ma->GetNE ());

    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    ScalarFE<eltyp, 1> faceint;

    // cout << "nbasis: " << nbasis << " ne: " << ma->GetNE() << " order: " <<
    // order << " number of ip: " << ir.Size() << endl;

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (ma);  // collection of tents in timeslab
    tps.PitchTents (dt, wavespeed); // adt = time slab height, wavespeed

    cout << "solving tents";
    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      // LocalHeap slh = lh.Split();  // split to threads
      HeapReset hr (lh);
      Tent *tent = tps.tents[tentnr];
      // cout << endl << "%%%% tent: " << i << " vert: " << tent->vertex << "
      // els: " << tent->els << endl;
      cout << *tent << endl;
      // if(tent->tbot==0 && tent->ttop-tent->tbot >= 0.19) cout <<
      // tent->ttop-tent->tbot << endl<< *tent << endl;

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;

      tel.SetCenter (center);
      tel.SetElSize (TentAdiam<D> (tent, ma));

      FlatMatrix<> elmat (nbasis, lh);
      FlatVector<> elvec (nbasis, lh);
      elmat = 0;
      elvec = 0;
      for (auto elnr : tent->els)
        {
          INT<D + 1> vnr = ma->GetEdgePNums (elnr);
          MappedIntegrationRule<D, D> mir (ir, ma->GetTrafo (elnr, lh),
                                           lh); // <dim  el, dim space>

          Mat<D + 1, D + 1> vtop = TentFaceVerts<D> (tent, elnr, ma, 1);
          Vec<D + 1> linearbasis_top = vtop.Row (D);
          Mat<D + 1, D + 1> vbot = TentFaceVerts<D> (tent, elnr, ma, 0);
          Vec<D + 1> linearbasis_bot = vbot.Row (D);
          for (int imip = 0; imip < mir.Size (); imip++)
            {
              Vec<D + 1> p;
              p.Range (0, D) = mir[imip].GetPoint ();

              // Integration over top of tent
              Vec<D + 1> n = TentFaceNormal<D> (vtop, 1);
              mir[imip].SetMeasure (TentFaceArea<D> (vtop));
              p (D) = faceint.Evaluate (ir[imip], linearbasis_top);

              FlatVector<> shape (nbasis, lh);
              FlatMatrix<> dshape (nbasis, D + 1, lh);

              tel.CalcShape (p, shape);
              tel.CalcDShape (p, dshape);

              for (int i = 0; i < nbasis; i++)
                {
                  for (int j = 0; j < nbasis; j++)
                    {
                      Vec<D> sig = -dshape.Row (i).Range (0, D);
                      Vec<D> tau = -dshape.Row (j).Range (0, D);
                      elmat (j, i)
                          += mir[imip].GetWeight ()
                             * (dshape (i, D) * dshape (j, D) * n (D)
                                    * (1 / (wavespeed * wavespeed))
                                + InnerProduct (sig, tau) * n (D)
                                + dshape (i, D)
                                      * InnerProduct (tau, n.Range (0, D))
                                + dshape (j, D)
                                      * InnerProduct (sig, n.Range (0, D)));
                    }
                }

              // Integration over bot of tent
              n = TentFaceNormal<D> (vbot, -1);
              mir[imip].SetMeasure (TentFaceArea<D> (vbot));
              p (D) = faceint.Evaluate (ir[imip], linearbasis_bot);

              tel.CalcShape (p, shape);
              tel.CalcDShape (p, dshape);

              int offset = elnr * ir.Size () * (D + 2) + imip * (D + 2);
              for (int j = 0; j < nbasis; j++)
                {
                  Vec<D> tau = -dshape.Row (j).Range (0, D);
                  elvec (j)
                      += mir[imip].GetWeight ()
                         * (-wavefront (offset + D + 1) * dshape (j, D) * n (D)
                                * (1 / (wavespeed * wavespeed))
                            - InnerProduct (
                                  wavefront.Range (offset + 1, offset + D + 1),
                                  tau)
                                  * n (D)
                            - wavefront (offset + D + 1)
                                  * InnerProduct (tau, n.Range (0, D))
                            - dshape (j, D)
                                  * InnerProduct (
                                      wavefront.Range (offset + 1,
                                                       offset + D + 1),
                                      n.Range (0, D))
                            + wavefront (offset) * shape (j));

                  for (int i = 0; i < nbasis; i++)
                    {
                      elmat (j, i)
                          += mir[imip].GetWeight () * (shape (i) * shape (j));
                    }
                }
            }
        } // close loop over tent elements

      // Integrate over side of tent
      for (auto surfel : ma->GetVertexSurfaceElements (tent->vertex))
        {
          auto sel_verts = ma->GetElVertices (ElementId (BND, surfel));
          Mat<D + 1, D + 1> v;
          v.Col (0).Range (0, D) = ma->GetPoint<D> (tent->vertex);
          v (D, 0) = tent->tbot;
          for (int n = 0; n < D; n++)
            {
              v.Col (n + 1).Range (0, D) = ma->GetPoint<D> (sel_verts[n]);
              v (D, n + 1) = tent->vertex == sel_verts[n]
                                 ? tent->ttop
                                 : tent->nbtime[tent->nbv.Pos (sel_verts[n])];
            }

          double A = TentFaceArea<D> (v);

          Vec<D + 1> n;
          n.Range (0, D)
              = TentFaceNormal<D - 1> (v.Cols (1, D).Rows (0, D - 1), 0);
          if (D == 1)
            n[0] = sgn_nozero<int> (tent->vertex - tent->nbv[0]);
          n (D) = 0;
          n[D] = 0;

          Mat<D + 1, D> map;
          for (int i = 0; i < D; i++)
            map.Col (i) = v.Col (i + 1) - v.Col (0);
          Vec<D + 1> shift = v.Col (0);

          for (int imip = 0; imip < ir.Size (); imip++)
            {
              Vec<D + 1> p = map * ir[imip].Point () + shift;
              FlatMatrix<> dshape (nbasis, D + 1, lh);
              tel.CalcDShape (p, dshape);
              p[D] += timeshift;
              double weight = A * ir[imip].Weight ();
              for (int j = 0; j < nbasis; j++)
                {
                  Vec<D> tau = -dshape.Row (j).Range (0, D);
                  elvec (j) -= weight * InnerProduct (tau, n.Range (0, D))
                               * TestSolution<D> (p, wavespeed)[D + 1];
                  for (int i = 0; i < nbasis; i++)
                    {
                      Vec<D> sig = -dshape.Row (i).Range (0, D);
                      elmat (j, i) += weight
                                      * InnerProduct (sig, n.Range (0, D))
                                      * dshape (j, D);
                    }
                }
            }
        }

      // solve
      Vector<> blavec = elvec;
      Matrix<> blamat = elmat;
      LapackSolve (elmat, elvec);
      Vector<> sol = elvec;
      if (L2Norm (blamat * sol - blavec) > 1e-7)
        cout << "error inverse: " << L2Norm (blamat * sol - blavec) << endl;

      // for(auto elnr: tent->els)
      // solutioncoeffs.Range(elnr*nbasis,(elnr+1)*nbasis) = sol;

      // eval solution on top of tent
      for (auto elnr : tent->els)
        {
          INT<D + 1> vnr = ma->GetEdgePNums (elnr);
          MappedIntegrationRule<D, D> mir (ir, ma->GetTrafo (elnr, lh),
                                           lh); // <dim  el, dim space>

          Mat<D + 1, D + 1> v = TentFaceVerts<D> (tent, elnr, ma, 1);
          Vec<D + 1> n = TentFaceNormal<D> (v, 1);
          Vec<D + 1> bs = v.Row (D);
          double A = TentFaceArea<D> (v);
          for (int imip = 0; imip < mir.Size (); imip++)
            {
              Vec<D + 1> p;
              p.Range (0, D) = mir[imip].GetPoint ();
              p (D) = faceint.Evaluate (ir[imip], bs);

              Matrix<> dshape (nbasis, D + 1);
              tel.CalcDShape (p, dshape);
              Vector<> shape (nbasis);
              tel.CalcShape (p, shape);

              int offset = elnr * ir.Size () * (D + 2) + imip * (D + 2);
              wavefront (offset) = InnerProduct (shape, sol);
              wavefront.Range (offset + 1, offset + D + 2)
                  = Trans (dshape) * sol;
              wavefront.Range (offset + 1, offset + D + 1) *= (-1);

              // if(imip==0 && L2Norm(wavefront.Range(offset,offset+D+2) -
              // TestSolution<D>(p,wavespeed))>pow(10,-order+2)) cout << " at "
              // << p << " error: " <<
              // L2Norm(wavefront.Range(offset,offset+D+2) -
              // TestSolution<D>(p,wavespeed)) << " isbndtent: " <<
              //ma->GetVertexSurfaceElements(tent->vertex).Size()<< " vert: "
              //<< tent->vertex <<  " tenth: " << tent->ttop-tent->tbot <<
              //endl;
            }
        }
    }); // end loop over tents
    cout << "...done" << endl;
  }

  // returns matrix where cols correspond to vertex coordinates of the
  // space-time element
  template <int D>
  Mat<D + 1, D + 1>
  TentFaceVerts (Tent *tent, int elnr, shared_ptr<MeshAccess> ma, bool top)
  {
    INT<D + 1> vnr = ma->GetElVertices (elnr);
    Mat<D + 1, D + 1> v;
    // determine linear basis function coeffs to use for tent face
    for (int ivert = 0; ivert < vnr.Size (); ivert++)
      {
        if (vnr[ivert] == tent->vertex)
          v (D, ivert) = top ? tent->ttop : tent->tbot;
        for (int k = 0; k < tent->nbv.Size (); k++)
          if (vnr[ivert] == tent->nbv[k])
            v (D, ivert) = tent->nbtime[k];
        v.Col (ivert).Range (0, D) = ma->GetPoint<D> (vnr[ivert]);
      }

    return v;
  }

  template <int D> double TentFaceArea (Mat<D + 1, D + 1> ve)
  {
    switch (D)
      {
      case 1:
        return L2Norm (ve.Col (0) - ve.Col (1));
        break;
      case 2:
        {
          double a = L2Norm (ve.Col (0) - ve.Col (1));
          double b = L2Norm (ve.Col (1) - ve.Col (2));
          double c = L2Norm (ve.Col (0) - ve.Col (2));
          SwapIfGreater<> (a, b);
          SwapIfGreater<> (a, c);
          SwapIfGreater<> (b, c);
          return 0.25
                 * sqrt ((a + (b + c)) * (c - (a - b)) * (c + (a - b))
                         * (a + (b - c)));
          break;
        }
      case 3:
        {
          double U = L2Norm (ve.Col (0) - ve.Col (1));
          double V = L2Norm (ve.Col (1) - ve.Col (2));
          double W = L2Norm (ve.Col (2) - ve.Col (0));
          double u = L2Norm (ve.Col (3) - ve.Col (2));
          double v = L2Norm (ve.Col (3) - ve.Col (0));
          double w = L2Norm (ve.Col (3) - ve.Col (1));

          double X = (w - U + v) * (U + v + w);
          double x = (U - v + w) * (v - w + U);
          double Y = (u - V + w) * (V + w + u);
          double y = (V - w + u) * (w - u + V);
          double Z = (v - W + u) * (W + u + v);
          double z = (W - u + v) * (u - v + W);

          double a = sqrt (x * Y * Z);
          double b = sqrt (y * Z * X);
          double c = sqrt (z * X * Y);
          double d = sqrt (x * y * z);

          return sqrt ((-a + b + c + d) * (a - b + c + d) * (a + b - c + d)
                       * (a + b + c - d))
                 / (192.0 * u * v * w);
        }
      }
  }

  template <int D> Vec<D + 1> TentFaceNormal (Mat<D + 1, D + 1> v, int top)
  {
    Vec<D + 1> normv;
    switch (D)
      {
      case 1:
        {
          normv (0) = v (1, 1) - v (1, 0);
          normv (1) = v (0, 0) - v (0, 1);
          break;
        }
      case 2:
        {
          Vec<D + 1> a = v.Col (0) - v.Col (1);
          Vec<D + 1> b = v.Col (0) - v.Col (2);
          normv (0) = a (1) * b (2) - a (2) * b (1);
          normv (1) = a (2) * b (0) - a (0) * b (2);
          normv (2) = a (0) * b (1) - a (1) * b (0);
          break;
        }
      case 3:
        {
          for (int d = 1; d <= D; d++)
            v.Col (d) = v.Col (0) - v.Col (d);

          for (unsigned int i = 0; i <= D; i++)
            {
              Mat<D, D> pS;
              for (unsigned int k = 0, c = 0; k < D + 1; k++)
                {
                  if (k == i)
                    continue;
                  pS.Row (c) = v.Row (k).Range (1, D + 1);
                  c++;
                }
              if ((i % 2) == 0)
                normv[i] = Det (pS);
              else
                normv[i] = -Det (pS);
            }
          break;
        }
      }
    normv /= L2Norm (normv);
    if (top == 1)
      normv *= sgn_nozero<double> (normv[D]);
    else if (top == -1)
      normv *= (-sgn_nozero<double> (normv[D]));
    return normv;
  }

  template <int D> Vec<D + 2> TestSolution (Vec<D + 1> p, double wavespeed)
  {
    double x = p[0];
    double t = p[D];
    Vec<D + 2> sol;
    int k = 3;
    if (D == 1)
      {
        sol[0] = sin (k * (wavespeed * t + x));
        sol[1] = -k * cos (k * (wavespeed * t + x));
        sol[2] = wavespeed * k * cos (k * (wavespeed * t + x));
      }
    else if (D == 2)
      {
        double y = p[1];
        double sq = sqrt (0.5);
        sol[0] = sin (wavespeed * t + sq * (x + y));
        sol[1] = -sq * cos (wavespeed * t + sq * (x + y));
        sol[2] = -sq * cos (wavespeed * t + sq * (x + y));
        sol[3] = wavespeed * cos (wavespeed * t + sq * (x + y));
      }
    else if (D == 3)
      {
        double y = p[1];
        double z = p[2];
        double sq = sqrt (1.0 / 3.0);
        sol[0] = sin (wavespeed * t + sq * (x + y + z));
        sol[1] = -sq * cos (wavespeed * t + sq * (x + y + z));
        sol[2] = -sq * cos (wavespeed * t + sq * (x + y + z));
        sol[3] = -sq * cos (wavespeed * t + sq * (x + y + z));
        sol[4] = wavespeed * cos (wavespeed * t + sq * (x + y + z));
      }
    return sol;
  }

  template <int D>
  Vector<> MakeWavefront (int order, shared_ptr<MeshAccess> ma,
                          double wavespeed, double time)
  {
    LocalHeap lh (10000000);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    Vector<> ic (ir.Size () * ma->GetNE () * (D + 2));
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        MappedIntegrationRule<D, D> mir (ir, ma->GetTrafo (elnr, lh),
                                         lh); // <dim  el, dim space>
        for (int imip = 0; imip < mir.Size (); imip++)
          {
            Vec<D + 1> p;
            p.Range (0, D) = mir[imip].GetPoint ();
            p[D] = time;
            int offset = elnr * ir.Size () * (D + 2) + imip * (D + 2);
            ic.Range (offset, offset + D + 2) = TestSolution<D> (p, wavespeed);
          }
      }
    return ic;
  }

  template <int D>
  double Postprocess (int order, shared_ptr<MeshAccess> ma, Vector<> wavefront,
                      Vector<> wavefront_corr)
  {
    double l2error = 0;
    LocalHeap lh (10000000);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        MappedIntegrationRule<D, D> mir (ir, ma->GetTrafo (elnr, lh),
                                         lh); // <dim  el, dim space>
        for (int imip = 0; imip < ir.Size (); imip++)
          {
            int offset = elnr * ir.Size () * (D + 2) + imip * (D + 2);
            l2error += (wavefront (offset) - wavefront_corr (offset))
                       * (wavefront (offset) - wavefront_corr (offset))
                       * mir[imip].GetWeight ();
          }
      }
    return sqrt (l2error);
  }

  template <typename T> void SwapIfGreater (T &a, T &b)
  {
    if (a < b)
      {
        T tmp (a);
        a = b;
        b = tmp;
      }
  }

  template <int D> double TentAdiam (Tent *tent, shared_ptr<MeshAccess> ma)
  {
    double anisotropicdiam = 0;
    int vnumber = tent->nbv.Size () + 2;

    Array<int> verts (vnumber);
    verts.Range (2, vnumber) = tent->nbv;
    verts[0] = tent->vertex;
    verts[1] = tent->vertex;

    Array<int> vtime (vnumber);
    vtime.Range (2, vnumber) = tent->nbtime;
    vtime[0] = tent->tbot;
    vtime[1] = tent->ttop;
    for (int k = 0; k < vnumber; k++)
      {
        for (int j = 0; j < vnumber; j++)
          {
            Vec<D> v1 = ma->GetPoint<D> (verts[j]);
            Vec<D> v2 = ma->GetPoint<D> (verts[k]);
            anisotropicdiam
                = max (anisotropicdiam, sqrt (L2Norm2 (v1 - v2)
                                              + pow (vtime[j] - vtime[k], 2)));
          }
      }
    return anisotropicdiam;
  }
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent (py::module m)
{
  m.def ("EvolveTents",
         [] (int order, shared_ptr<MeshAccess> ma, double wavespeed, double dt,
             Vector<> wavefront,
             double timeshift = 0) -> Vector<> //-> shared_ptr<MeshAccess>
         {
           int D = ma->GetDimension ();
           if (D == 1)
             EvolveTents<1> (order, ma, wavespeed, dt, wavefront, timeshift);
           else if (D == 2)
             EvolveTents<2> (order, ma, wavespeed, dt, wavefront, timeshift);
           else if (D == 3)
             EvolveTents<3> (order, ma, wavespeed, dt, wavefront, timeshift);
           return wavefront;
         } //, py::call_guard<py::gil_scoped_release>()
  );
  m.def ("EvolveTentsMakeWavefront",
         [] (int order, shared_ptr<MeshAccess> ma, double wavespeed,
             double time) -> Vector<> //-> shared_ptr<MeshAccess>
         {
           int D = ma->GetDimension ();
           Vector<> wavefront;
           if (D == 1)
             wavefront = MakeWavefront<1> (order, ma, wavespeed, time);
           else if (D == 2)
             wavefront = MakeWavefront<2> (order, ma, wavespeed, time);
           else if (D == 3)
             wavefront = MakeWavefront<3> (order, ma, wavespeed, time);
           return wavefront;
         });
  m.def ("EvolveTentsPostProcess",
         [] (int order, shared_ptr<MeshAccess> ma, Vector<> wavefront,
             Vector<> wavefront_corr) -> double {
           int D = ma->GetDimension ();
           double l2error;
           if (D == 1)
             l2error = Postprocess<1> (order, ma, wavefront, wavefront_corr);
           else if (D == 2)
             l2error = Postprocess<2> (order, ma, wavefront, wavefront_corr);
           else if (D == 3)
             l2error = Postprocess<3> (order, ma, wavefront, wavefront_corr);
           return l2error;
         });
}
#endif // NGS_PYTHON
