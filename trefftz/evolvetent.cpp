#include "evolvetent.hpp"
#include "trefftzwavefe.hpp"
#include "tents/tents.hpp"
#include "testcases.hpp"
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
                    double dt, SliceMatrix<double> wavefront, double timeshift)
  {
    LocalHeap lh (100000000);

    int nsimd = SIMD<double>::Size ();

    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    int nip = ir.Size ();

    SIMD_IntegrationRule sir (eltyp, order * 2);
    int snip = sir.Size () * nsimd;

    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces
    TrefftzWaveFE<D + 1> tel (order, wavespeed);
    int nbasis = tel.GetNBasis ();

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (ma);      // collection of tents in timeslab
    tps.PitchTents (dt, wavespeed + 1); // adt = time slab height, wavespeed

    // FlatVector<SIMD<double>> simd_wavefront(sir.Size()*ma->GetNE()*(D+2),
    // &wavefront[0]);

    cout << "solving tents";
    static Timer ttent ("tent", 2);
    static Timer tint ("tentint", 2);
    static Timer tcalcshape ("tentcalcshape", 2);
    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      HeapReset hr (lh);
      Tent *tent = tps.tents[tentnr];

      RegionTimer reg (ttent);

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      tel.SetCenter (center);
      tel.SetElSize (TentAdiam<D> (tent, ma));

      FlatMatrix<> elmat (nbasis, lh);
      FlatVector<> elvec (nbasis, lh);
      elmat = 0;
      elvec = 0;

      LocalHeap slh = lh.Split (); // split to threads
      for (auto elnr : tent->els)
        {
          HeapReset hr (slh);

          SIMD_MappedIntegrationRule<D, D + 1> smir (
              sir, ma->GetTrafo (elnr, slh), slh);
          SIMD_MappedIntegrationRule<D, D> smir_fix (
              sir, ma->GetTrafo (elnr, slh), slh);
          for (int imip = 0; imip < sir.Size (); imip++)
            smir[imip].Point ().Range (0, D)
                = smir_fix[imip].Point ().Range (0, D);

          Mat<D + 1> vert;     // vertices of tent face
          Vec<D + 1> linbasis; // coeffs for linear face fct
          Mat<D + 1> Dmat;
          FlatVector<SIMD<double>> mirtimes (sir.Size (), slh);

          /// Integration over top of tent
          vert = TentFaceVerts<D> (tent, elnr, ma, 1);
          linbasis = vert.Row (D);
          faceint.Evaluate (sir, linbasis, mirtimes);
          for (int imip = 0; imip < sir.Size (); imip++)
            smir[imip].Point () (D) = mirtimes[imip];

          FlatMatrix<SIMD<double>> simddshapes ((D + 1) * nbasis, sir.Size (),
                                                slh);
          tel.CalcDShape (smir, simddshapes);
          FlatMatrix<> bbmat (nbasis, (D + 1) * snip, &simddshapes (0, 0)[0]);

          TentDmat<D> (Dmat, vert, 1, wavespeed);
          FlatMatrix<> bdbmat ((D + 1) * snip, nbasis, slh);
          bdbmat = 0;
          for (int imip = 0; imip < snip; imip++)
            for (int r = 0; r < (D + 1); r++)
              for (int d = 0; d < D + 1; d++)
                bdbmat.Row (r * snip + imip)
                    += Dmat (r, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                       * bbmat.Col (d * snip + imip);
          elmat += bbmat * bdbmat;

          /// Integration over bot of tent
          vert = TentFaceVerts<D> (tent, elnr, ma, 0);
          linbasis = vert.Row (D);
          faceint.Evaluate (sir, linbasis, mirtimes);
          for (int imip = 0; imip < sir.Size (); imip++)
            smir[imip].Point () (D) = mirtimes[imip];

          FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
          tel.CalcShape (smir, simdshapes);
          tel.CalcDShape (smir, simddshapes);

          TentDmat<D> (Dmat, vert, -1, wavespeed);
          FlatVector<> bdbvec ((D + 1) * snip, slh);
          bdbvec = 0;
          for (int imip = 0; imip < snip; imip++)
            for (int r = 0; r < (D + 1); r++)
              for (int d = 0; d < D + 1; d++)
                bdbvec (r * snip + imip)
                    += Dmat (r, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                       * wavefront (elnr, nip + (imip % nip) * (D + 1) + d);
          elvec -= bbmat * bdbvec;

          // stabilization to recover second order solution
          for (int imip = 0; imip < sir.Size (); imip++)
            simdshapes.Col (imip)
                *= sqrt (TentFaceArea<D> (vert)) * sqrt (sir[imip].Weight ());
          AddABt (simdshapes, simdshapes, elmat);
          for (int imip = 0; imip < sir.Size (); imip++)
            simdshapes.Col (imip)
                *= sqrt (TentFaceArea<D> (vert)) * sqrt (sir[imip].Weight ());
          FlatMatrix<> shapes (nbasis, sir.Size () * nsimd,
                               &simdshapes (0, 0)[0]);
          elvec += shapes * wavefront.Row (elnr).Range (0, nip);
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

          Vec<D + 1> n;
          n.Range (0, D)
              = TentFaceNormal<D> (v.Cols (1, D + 1).Rows (0, D), 0);
          if (D == 1)
            n[0] = sgn_nozero<int> (tent->vertex - tent->nbv[0]);
          n[D] = 0;

          Mat<D + 1, D> map;
          for (int i = 0; i < D; i++)
            map.Col (i) = v.Col (i + 1) - v.Col (0);
          Vec<D + 1> shift = v.Col (0);

          SIMD_MappedIntegrationRule<D, D + 1> smir (
              sir, ma->GetTrafo (0, slh), slh);
          for (int imip = 0; imip < ir.Size (); imip++)
            smir[imip].Point ()
                = map * sir[imip].operator Vec<D, SIMD<double>> () + shift;

          FlatMatrix<SIMD<double>> simddshapes ((D + 1) * nbasis, sir.Size (),
                                                slh);
          tel.CalcDShape (smir, simddshapes);
          FlatMatrix<double> bbmat (nbasis, (D + 1) * snip,
                                    &simddshapes (0, 0)[0]);

          Mat<D + 1> Dmat = 0;
          Dmat.Row (D).Range (0, D) = -TentFaceArea<D> (v) * n.Range (0, D);
          FlatMatrix<double> bdbmat ((D + 1) * snip, nbasis, slh);
          bdbmat = 0;
          for (int imip = 0; imip < snip; imip++)
            for (int r = 0; r < (D + 1); r++)
              for (int d = 0; d < D + 1; d++)
                bdbmat.Row (r * snip + imip)
                    += Dmat (r, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                       * bbmat.Col (d * snip + imip);

          elmat += bbmat * bdbmat;

          Vector<> bc = EvalBC<D> (smir, wavespeed, timeshift);
          FlatVector<> bdbvec ((D + 1) * snip, slh);
          bdbvec = 0;
          for (int imip = 0; imip < snip; imip++)
            for (int r = 0; r < (D + 1); r++)
              for (int d = 0; d < D + 1; d++)
                // use Dmat transposed
                bdbvec (r * snip + imip)
                    += Dmat (d, r) * sir[imip / nsimd].Weight ()[imip % nsimd]
                       * bc ((imip % nip) * (D + 1) + d);

          elvec -= bbmat * bdbvec;
        }

      // solve
      LapackSolve (elmat, elvec);
      FlatVector<> sol (nbasis, &elvec (0));

      double tenterror = 0;
      // eval solution on top of tent
      for (auto elnr : tent->els)
        {
          MappedIntegrationRule<D, D> mir (ir, ma->GetTrafo (elnr, slh),
                                           slh); // <dim  el, dim space>

          Mat<D + 1, D + 1> v = TentFaceVerts<D> (tent, elnr, ma, 1);
          Vec<D + 1> bs = v.Row (D);
          for (int imip = 0; imip < nip; imip++)
            mir[imip].Point () (D) = faceint.Evaluate (ir[imip], bs);

          FlatMatrix<> shapes (nbasis, nip, slh);
          FlatMatrix<> dshapes (nbasis, (D + 1) * nip, slh);
          tel.CalcDShape (mir, dshapes);
          tel.CalcShape (mir, shapes);

          wavefront.Row (elnr).Range (0, nip) = Trans (shapes) * sol;
          wavefront.Row (elnr).Range (nip, nip + nip * (D + 1))
              = Trans (dshapes) * sol;
          // p[D] += timeshift;
          // tenterror +=
          // (wavefront(offset)-TestSolution<D>(p,wavespeed)[0])*(wavefront(offset)-TestSolution<D>(p,wavespeed)[0])*ir[imip].Weight()
          // * A;
        }
      // cout << "error tent: " << sqrt(tenterror) << endl;
    }); // end loop over tents
    cout << "...done" << endl;
  }

  template <int D>
  void TentDmat (Mat<D + 1> &Dmat, Mat<D + 1> v, int top, double wavespeed)
  {
    Vec<D + 1> n = TentFaceNormal<D + 1> (v, top);
    Dmat = n (D) * Id<D + 1> ();
    Dmat.Row (D).Range (0, D) = -n.Range (0, D);
    Dmat.Col (D).Range (0, D) = -n.Range (0, D);
    Dmat (D, D) *= 1.0 / (wavespeed * wavespeed);
    Dmat *= TentFaceArea<D> (v);
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
        else
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

  template <int D> Vec<D> TentFaceNormal (Mat<D, D> v, int top)
  {
    Vec<D> normv;
    switch (D)
      {
      case 2:
        {
          normv (0) = v (1, 1) - v (1, 0);
          normv (1) = v (0, 0) - v (0, 1);
          break;
        }
      case 3:
        {
          Vec<D + 1> a = v.Col (0) - v.Col (1);
          Vec<D + 1> b = v.Col (0) - v.Col (2);
          normv (0) = a (1) * b (2) - a (2) * b (1);
          normv (1) = a (2) * b (0) - a (0) * b (2);
          normv (2) = a (0) * b (1) - a (1) * b (0);
          break;
        }
      case 4:
        {
          for (int d = 1; d < D; d++)
            v.Col (d) = v.Col (0) - v.Col (d);

          for (unsigned int i = 0; i < D; i++)
            {
              Mat<D - 1, D - 1> pS;
              for (unsigned int k = 0, c = 0; k < D; k++)
                {
                  if (k == i)
                    continue;
                  pS.Row (c) = v.Row (k).Range (1, D);
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
      normv *= sgn_nozero<double> (normv[D - 1]);
    else if (top == -1)
      normv *= (-sgn_nozero<double> (normv[D - 1]));
    return normv;
  }

  template <int D>
  Vector<> EvalBC (const SIMD_MappedIntegrationRule<D, D + 1> &mir,
                   double wavespeed, double timeshift)
  {
    int nsimd = SIMD<double>::Size ();
    Vector<> bc ((D + 1) * mir.Size () * nsimd);
    for (int imip = 0; imip < mir.Size (); imip++)
      {
        Vector<SIMD<double>> sp = mir[imip].GetPoint ();
        sp[D] += timeshift;
        Vec<D + 1> p;
        for (int s = 0; s < nsimd; s++)
          {
            for (int d = 0; d < D + 1; d++)
              p[d] = sp[d][s];
            bc.Range (imip * (D + 1) * nsimd + s * (D + 1),
                      (imip) * (D + 1) * nsimd + (s + 1) * (D + 1))
                = TestSolution<D> (p, wavespeed).Range (1, D + 2);
          }
      }
    return bc;
  }

  template <int D>
  Matrix<> MakeWavefront (int order, shared_ptr<MeshAccess> ma,
                          double wavespeed, double time)
  {
    LocalHeap lh (10000000);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    int nip = ir.Size ();
    Matrix<> ic (ma->GetNE (), nip * (D + 2));
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        MappedIntegrationRule<D, D + 1> mir (ir, ma->GetTrafo (elnr, lh),
                                             lh); // <dim  el, dim space>
        for (int imip = 0; imip < nip; imip++)
          {
            mir[imip].Point () (D) = time;
            ic (elnr, imip)
                = TestSolution<D> (mir[imip].Point (), wavespeed)[0];
            ic.Row (elnr).Range (nip + imip * (D + 1),
                                 nip + (imip + 1) * (D + 1))
                = TestSolution<D> (mir[imip].Point (), wavespeed)
                      .Range (1, D + 2);
          }
      }
    return ic;
  }

  template <int D>
  double Postprocess (int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront,
                      Matrix<> wavefront_corr)
  {
    double l2error = 0;
    LocalHeap lh (10000000);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        for (int imip = 0; imip < ir.Size (); imip++)
          {
            l2error += (wavefront (elnr, imip) - wavefront_corr (elnr, imip))
                       * (wavefront (elnr, imip) - wavefront_corr (elnr, imip))
                       * ir[imip].Weight ();
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
             Matrix<> wavefront,
             double timeshift = 0) -> Matrix<> //-> shared_ptr<MeshAccess>
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
             double time) -> Matrix<> //-> shared_ptr<MeshAccess>
         {
           int D = ma->GetDimension ();
           Matrix<> wavefront;
           if (D == 1)
             wavefront = MakeWavefront<1> (order, ma, wavespeed, time);
           else if (D == 2)
             wavefront = MakeWavefront<2> (order, ma, wavespeed, time);
           else if (D == 3)
             wavefront = MakeWavefront<3> (order, ma, wavespeed, time);
           return wavefront;
         });
  m.def ("EvolveTentsPostProcess",
         [] (int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront,
             Matrix<> wavefront_corr) -> double {
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
