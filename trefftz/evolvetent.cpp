#include "evolvetent.hpp"

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
                    double dt, SliceMatrix<double> wavefront, double timeshift,
                    shared_ptr<CoefficientFunction> bddatum)
  {
    LocalHeap lh (100000000);

    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    int nsimd = SIMD<double>::Size ();
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int snip = sir.Size () * nsimd;

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (ma);      // collection of tents in timeslab
    tps.PitchTents (dt, wavespeed + 1); // adt = time slab height, wavespeed

    cout << "solving " << tps.tents.Size () << " tents ";
    static Timer ttent ("tent", 2);
    static Timer ttentel ("tentel", 2);
    static Timer ttentbnd ("tentbnd", 2);
    static Timer ttenteval ("tenteval", 2);

    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      RegionTimer reg (ttent);
      LocalHeap slh = lh.Split (); // split to threads
      HeapReset hr (slh);
      Tent *tent = tps.tents[tentnr];

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      TrefftzWaveFE<D + 1> tel (order, wavespeed, center,
                                TentAdiam<D> (tent, ma, wavespeed));
      int nbasis = tel.GetNBasis ();

      FlatMatrix<> elmat (nbasis, slh);
      FlatVector<> elvec (nbasis, slh);
      elmat = 0;
      elvec = 0;

      FlatMatrix<SIMD<double>> *topdshapes[tent->els.Size ()];
      for (auto &tds : topdshapes)
        tds = new FlatMatrix<SIMD<double>> ((D + 1) * nbasis, sir.Size (),
                                            slh);

      // Integrate top and bottom space-like tent faces
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          RegionTimer reg1 (ttentel);
          CalcTentEl<D> (tent->els[elnr], tent, tel, ma, wavefront, sir, slh,
                         elmat, elvec, *topdshapes[elnr]);
        }

      // Integrate boundary tent
      for (auto surfel : ma->GetVertexSurfaceElements (tent->vertex))
        {
          RegionTimer reg2 (ttentbnd);
          CalcTentBndEl<D> (surfel, tent, tel, ma, bddatum, timeshift, sir,
                            slh, elmat, elvec);
        }

      // solve
      LapackSolve (elmat, elvec);
      FlatVector<> sol (nbasis, &elvec (0));

      // modifications for first order solver
      // Matrix<> elmat2 = elmat.Rows(1,nbasis).Cols(1,nbasis);
      // Vector<> elvec2 = elvec.Range(1,nbasis);
      // LapackSolve(elmat2,elvec2);
      // Vector<> sol(nbasis);
      // sol[0] = 0;
      // for(int i = 1; i<nbasis;i++)
      // sol[i]=elvec2[i-1];

      // eval solution on top of tent
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          RegionTimer reg3 (ttenteval);
          CalcTentElEval<D> (tent->els[elnr], tent, tel, ma, wavefront, sir,
                             slh, sol, *topdshapes[elnr]);
        }
    }); // end loop over tents
    cout << "...done" << endl;
  }

  template <int D>
  void EvolveTents (int order, shared_ptr<MeshAccess> ma, Vector<> wavespeed,
                    double dt, SliceMatrix<double> wavefront, double timeshift,
                    shared_ptr<CoefficientFunction> bddatum)
  {
    LocalHeap lh (100000000);

    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    int nsimd = SIMD<double>::Size ();
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int snip = sir.Size () * nsimd;
    int ndomains = ma->GetNDomains ();

    double max_wavespeed = wavespeed[0];
    for (double c : wavespeed)
      max_wavespeed = max (c, max_wavespeed);

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (ma); // collection of tents in timeslab
    tps.PitchTents (dt,
                    max_wavespeed + 1); // adt = time slab height, wavespeed

    cout << "solving " << tps.tents.Size () << " tents ";
    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      LocalHeap slh = lh.Split (); // split to threads
      HeapReset hr (slh);
      Tent *tent = tps.tents[tentnr];

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      TrefftzWaveFE<D + 1> tel (order, max_wavespeed, center,
                                TentAdiam<D> (tent, ma, wavespeed[0]));
      int nbasis = tel.GetNBasis ();

      // check if tent vertex is on boundary between domains
      int ndomains = 1;
      for (auto surfel : ma->GetVertexSurfaceElements (tent->vertex))
        {
          int in;
          int out;
          ma->GetSElNeighbouringDomains (surfel, in, out);
          if (out != 0)
            ndomains = 2;
        }

      FlatMatrix<> elmat (ndomains * nbasis, slh);
      FlatVector<> elvec (ndomains * nbasis, slh);
      elmat = 0;
      elvec = 0;

      FlatMatrix<SIMD<double>> *topdshapes[tent->els.Size ()];
      for (auto &tds : topdshapes)
        tds = new FlatMatrix<SIMD<double>> ((D + 1) * nbasis, sir.Size (),
                                            slh);

      // Integrate top and bottom space-like tent faces
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          int eli = ma->GetElIndex (ElementId (VOL, tent->els[elnr]));
          tel.SetWavespeed (wavespeed[eli]);
          if (ndomains == 1) // tent vertex is inside a domain
            eli = 0;
          SliceMatrix<> subm = elmat.Cols (eli * nbasis, (eli + 1) * nbasis)
                                   .Rows (eli * nbasis, (eli + 1) * nbasis);
          SliceVector<> subv = elvec.Range (eli * nbasis, (eli + 1) * nbasis);
          CalcTentEl<D> (tent->els[elnr], tent, tel, ma, wavefront, sir, slh,
                         subm, subv, *topdshapes[elnr]);
        }

      // Integrate boundary tent
      for (auto surfel : ma->GetVertexSurfaceElements (tent->vertex))
        {
          int in;
          int out;
          ma->GetSElNeighbouringDomains (surfel, in, out);

          if (out == 0)
            {
              tel.SetWavespeed (wavespeed[in - 1]);
              if (ndomains == 1) // tent vertex is inside a domain
                in = 1;
              SliceMatrix<> subm = elmat.Cols ((in - 1) * nbasis, in * nbasis)
                                       .Rows ((in - 1) * nbasis, in * nbasis);
              SliceVector<> subv
                  = elvec.Range ((in - 1) * nbasis, in * nbasis);
              CalcTentBndEl<D> (surfel, tent, tel, ma, bddatum, timeshift, sir,
                                slh, subm, subv);
            }
          else
            {
              int nbasis = tel.GetNBasis ();
              int nsimd = SIMD<double>::Size ();
              int snip = sir.Size () * nsimd;

              // get vertices of tent face
              Mat<D + 1> vert = TentFaceVerts<D> (tent, surfel, ma, 0);
              // build normal vector
              Vec<D + 1> n;
              n.Range (0, D)
                  = TentFaceNormal<D> (vert.Cols (1, D + 1).Rows (0, D), 0);
              if (D == 1) // D=1 special case
                n[0] = sgn_nozero<int> (tent->vertex - tent->nbv[0]);
              n[D] = 0; // time-like faces only
              // build mapping to physical boundary simplex
              Mat<D + 1, D> map;
              for (int i = 0; i < D; i++)
                map.Col (i) = vert.Col (i + 1) - vert.Col (0);
              Vec<D + 1> shift = vert.Col (0);

              SIMD_MappedIntegrationRule<D, D + 1> smir (
                  sir, ma->GetTrafo (0, slh), -1, slh);
              for (int imip = 0; imip < sir.Size (); imip++)
                smir[imip].Point ()
                    = map * sir[imip].operator Vec<D, SIMD<double>> () + shift;

              tel.SetWavespeed (wavespeed[in - 1]);
              FlatMatrix<SIMD<double>> simddshapes1 ((D + 1) * nbasis,
                                                     sir.Size (), slh);
              tel.CalcDShape (smir, simddshapes1);
              FlatMatrix<double> bbmat1 (nbasis, (D + 1) * snip,
                                         &simddshapes1 (0, 0)[0]);

              tel.SetWavespeed (wavespeed[out - 1]);
              FlatMatrix<SIMD<double>> simddshapes2 ((D + 1) * nbasis,
                                                     sir.Size (), slh);
              tel.CalcDShape (smir, simddshapes2);
              FlatMatrix<double> bbmat2 (nbasis, (D + 1) * snip,
                                         &simddshapes2 (0, 0)[0]);

              Mat<D + 1> Dmat1 = 0;
              Dmat1.Row (D).Range (0, D)
                  = -0.5 * TentFaceArea<D> (vert)
                    * n.Range (0, D); // 0.5 for DG average
              Dmat1.Col (D).Range (0, D)
                  = -0.5 * TentFaceArea<D> (vert) * n.Range (0, D);
              Mat<D + 1> Dmat2 = -1 * Dmat1;

              FlatMatrix<> *bdbmat[4];
              for (int i = 0; i < 4; i++)
                {
                  bdbmat[i] = new FlatMatrix<> ((D + 1) * snip, nbasis, slh);
                  *bdbmat[i] = 0;
                }
              // double alpha = 0.5;

              for (int imip = 0; imip < snip; imip++)
                for (int r = 0; r < (D + 1); r++)
                  for (int d = 0; d < D + 1; d++)
                    {
                      bdbmat[0]->Row (r * snip + imip)
                          += Dmat1 (r, d)
                             * sir[imip / nsimd].Weight ()[imip % nsimd]
                             * bbmat1.Col (d * snip + imip);
                      bdbmat[1]->Row (r * snip + imip)
                          += Dmat1 (r, d)
                             * sir[imip / nsimd].Weight ()[imip % nsimd]
                             * bbmat2.Col (d * snip + imip);
                      bdbmat[2]->Row (r * snip + imip)
                          += Dmat2 (r, d)
                             * sir[imip / nsimd].Weight ()[imip % nsimd]
                             * bbmat1.Col (d * snip + imip);
                      bdbmat[3]->Row (r * snip + imip)
                          += Dmat2 (r, d)
                             * sir[imip / nsimd].Weight ()[imip % nsimd]
                             * bbmat2.Col (d * snip + imip);
                    }
              elmat.Cols ((in - 1) * nbasis, in * nbasis)
                  .Rows ((in - 1) * nbasis, in * nbasis)
                  += bbmat1 * (*bdbmat[0]);
              elmat.Cols ((out - 1) * nbasis, out * nbasis)
                  .Rows ((in - 1) * nbasis, in * nbasis)
                  += bbmat1 * (*bdbmat[1]);
              elmat.Cols ((in - 1) * nbasis, in * nbasis)
                  .Rows ((out - 1) * nbasis, out * nbasis)
                  += bbmat2 * (*bdbmat[2]);
              elmat.Cols ((out - 1) * nbasis, out * nbasis)
                  .Rows ((out - 1) * nbasis, out * nbasis)
                  += bbmat2 * (*bdbmat[3]);
            }
        }

      // solve
      LapackSolve (elmat, elvec);
      FlatVector<> sol (ndomains * nbasis, &elvec (0));

      // eval solution on top of tent
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          int eli = ma->GetElIndex (ElementId (VOL, tent->els[elnr]));
          tel.SetWavespeed (wavespeed[eli]);
          if (ndomains == 1)
            eli = 0;
          CalcTentElEval<D> (tent->els[elnr], tent, tel, ma, wavefront, sir,
                             slh, sol.Range (eli * nbasis, (eli + 1) * nbasis),
                             *topdshapes[elnr]);
        }
    }); // end loop over tents
    cout << "...done" << endl;
  }

  template <int D>
  void
  CalcTentEl (int elnr, Tent *tent, TrefftzWaveFE<D + 1> tel,
              shared_ptr<MeshAccess> ma, SliceMatrix<double> &wavefront,
              SIMD_IntegrationRule &sir, LocalHeap &slh, SliceMatrix<> elmat,
              SliceVector<> elvec, SliceMatrix<SIMD<double>> simddshapes)
  {
    static Timer tint1 ("tent int calcshape", 2);
    static Timer tint2 ("tent int mat&vec", 2);

    HeapReset hr (slh);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    double wavespeed = tel.GetWavespeed ();
    int nbasis = tel.GetNBasis ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces

    SIMD_MappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (elnr, slh),
                                               -1, slh);
    SIMD_MappedIntegrationRule<D, D> smir_fix (sir, ma->GetTrafo (elnr, slh),
                                               slh);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ().Range (0, D) = smir_fix[imip].Point ().Range (0, D);

    Mat<D + 1> vert;     // vertices of tent face
    Vec<D + 1> linbasis; // coeffs for linear face fct
    Mat<D + 1> Dmat;
    FlatVector<SIMD<double>> mirtimes (sir.Size (), slh);

    /// Integration over bot of tent
    vert = TentFaceVerts<D> (tent, elnr, ma, -1);
    linbasis = vert.Row (D);
    faceint.Evaluate (sir, linbasis, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    tint1.Start ();
    FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
    tel.CalcShape (smir, simdshapes);
    tel.CalcDShape (smir, simddshapes);
    FlatMatrix<> bbmat (nbasis, (D + 1) * snip, &simddshapes (0, 0)[0]);
    tint1.Stop ();

    tint2.Start ();
    TentDmat<D> (Dmat, vert, -1, wavespeed);
    FlatVector<> bdbvec ((D + 1) * snip, slh);
    bdbvec = 0;
    for (int imip = 0; imip < snip; imip++)
      for (int r = 0; r < (D + 1); r++)
        for (int d = 0; d < D + 1; d++)
          bdbvec (r * snip + imip)
              += Dmat (r, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                 * wavefront (elnr, snip + d * snip + imip);
    elvec -= bbmat * bdbvec;

    // stabilization to recover second order solution
    for (int imip = 0; imip < sir.Size (); imip++)
      simdshapes.Col (imip)
          *= sqrt (TentFaceArea<D> (vert) * sir[imip].Weight ());
    AddABt (simdshapes, simdshapes, elmat);
    for (int imip = 0; imip < sir.Size (); imip++)
      simdshapes.Col (imip)
          *= sqrt (TentFaceArea<D> (vert) * sir[imip].Weight ());
    FlatMatrix<> shapes (nbasis, snip, &simdshapes (0, 0)[0]);
    elvec += shapes * wavefront.Row (elnr).Range (0, snip);
    tint2.Stop ();

    /// Integration over top of tent
    vert = TentFaceVerts<D> (tent, elnr, ma, 1);
    linbasis = vert.Row (D);
    faceint.Evaluate (sir, linbasis, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    tint1.Start ();
    // FlatMatrix<SIMD<double>> simddshapes((D+1)*nbasis,sir.Size(),slh);
    tel.CalcDShape (smir, simddshapes);
    tint1.Stop ();

    tint2.Start ();
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
    tint2.Stop ();
  }

  template <int D>
  void CalcTentBndEl (int surfel, Tent *tent, TrefftzWaveFE<D + 1> tel,
                      shared_ptr<MeshAccess> ma,
                      shared_ptr<CoefficientFunction> bddatum,
                      double timeshift, SIMD_IntegrationRule &sir,
                      LocalHeap &slh, SliceMatrix<> elmat, SliceVector<> elvec)
  {
    HeapReset hr (slh);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    double wavespeed = tel.GetWavespeed ();
    int nbasis = tel.GetNBasis ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;

    // get vertices of tent face
    Mat<D + 1> vert = TentFaceVerts<D> (tent, surfel, ma, 0);
    // build normal vector
    Vec<D + 1> n;
    n.Range (0, D) = TentFaceNormal<D> (vert.Cols (1, D + 1).Rows (0, D), 0);
    if (D == 1) // D=1 special case
      n[0] = sgn_nozero<int> (tent->vertex - tent->nbv[0]);
    n[D] = 0; // time-like faces only
    // build mapping to physical boundary simplex
    Mat<D + 1, D> map;
    for (int i = 0; i < D; i++)
      map.Col (i) = vert.Col (i + 1) - vert.Col (0);
    Vec<D + 1> shift = vert.Col (0);

    SIMD_MappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (0, slh), -1,
                                               slh);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ()
          = map * sir[imip].operator Vec<D, SIMD<double>> () + shift;

    FlatMatrix<SIMD<double>> simddshapes ((D + 1) * nbasis, sir.Size (), slh);
    tel.CalcDShape (smir, simddshapes);
    FlatMatrix<double> bbmat (nbasis, (D + 1) * snip, &simddshapes (0, 0)[0]);

    Mat<D + 1> Dmat = 0;
    Dmat.Row (D).Range (0, D) = -TentFaceArea<D> (vert) * n.Range (0, D);
    FlatMatrix<double> bdbmat ((D + 1) * snip, nbasis, slh);
    bdbmat = 0;

    for (int imip = 0; imip < smir.Size (); imip++)
      smir[imip].Point ()[D] += timeshift;
    FlatMatrix<SIMD<double>> bdeval ((D + 2), sir.Size (), slh);
    bdeval = 0;
    bddatum->Evaluate (smir, bdeval);

    FlatVector<> bdbvec ((D + 1) * snip, slh);
    bdbvec = 0;
    if (D > 1 && ma->GetMaterial (ElementId (BND, surfel)) == "neumann")
      {
        for (int imip = 0; imip < snip; imip++)
          for (int r = 0; r < (D + 1); r++)
            bdbmat.Row (r * snip + imip)
                += Dmat (D, r) * sir[imip / nsimd].Weight ()[imip % nsimd]
                   * bbmat.Col (D * snip + imip); // neumann
        elmat += bbmat * bdbmat;

        for (int imip = 0; imip < snip; imip++)
          for (int d = 0; d < D + 1; d++)
            bdbvec (D * snip + imip)
                += Dmat (D, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                   * bdeval (d + 1, imip / nsimd)[imip % nsimd];
        elvec -= bbmat * bdbvec;
      }
    else
      { // dirichlet
        double alpha = 0.5;
        Dmat (D, D) = TentFaceArea<D> (vert) * alpha;
        for (int imip = 0; imip < snip; imip++)
          // for(int r=0;r<(D+1);r++) // r=D since only last row non-zero
          for (int d = 0; d < D + 1; d++)
            bdbmat.Row (D * snip + imip)
                += Dmat (D, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                   * bbmat.Col (d * snip + imip);
        elmat += bbmat * bdbmat;

        Dmat (D, D) *= -1.0;
        for (int imip = 0; imip < snip; imip++)
          for (int r = 0; r < (D + 1); r++)
            bdbvec (r * snip + imip)
                += Dmat (D, r) * sir[imip / nsimd].Weight ()[imip % nsimd]
                   * bdeval (D + 1,
                             imip
                                 / nsimd)[imip % nsimd]; // use Dmat transposed
        elvec -= bbmat * bdbvec;
      }
  }

  template <int D>
  void
  CalcTentElEval (int elnr, Tent *tent, TrefftzWaveFE<D + 1> tel,
                  shared_ptr<MeshAccess> ma, SliceMatrix<> &wavefront,
                  SIMD_IntegrationRule &sir, LocalHeap &slh, SliceVector<> sol,
                  SliceMatrix<SIMD<double>> simddshapes)
  {
    HeapReset hr (slh);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    int nbasis = tel.GetNBasis ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces

    SIMD_MappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (elnr, slh),
                                               -1, slh);
    SIMD_MappedIntegrationRule<D, D> smir_fix (sir, ma->GetTrafo (elnr, slh),
                                               slh);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ().Range (0, D) = smir_fix[imip].Point ().Range (0, D);

    Mat<D + 1> v = TentFaceVerts<D> (tent, elnr, ma, 1);
    Vec<D + 1> bs = v.Row (D);
    FlatVector<SIMD<double>> mirtimes (sir.Size (), slh);
    faceint.Evaluate (sir, bs, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
    // FlatMatrix<SIMD<double>> simddshapes((D+1)*nbasis,sir.Size(),slh);
    tel.CalcShape (smir, simdshapes);
    // tel.CalcDShape(smir,simddshapes);
    FlatMatrix<> dshapes (nbasis, (D + 1) * snip, &simddshapes (0, 0)[0]);
    FlatMatrix<> shapes (nbasis, snip, &simdshapes (0, 0)[0]);

    wavefront.Row (elnr).Range (0, snip) = Trans (shapes.Cols (0, snip)) * sol;
    wavefront.Row (elnr).Range (snip, snip + snip * (D + 1))
        = Trans (dshapes) * sol;
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
  TentFaceVerts (Tent *tent, int elnr, shared_ptr<MeshAccess> ma, int top)
  {
    Mat<D + 1, D + 1> v;
    if (top == 0) // boundary element
      {
        auto sel_verts = ma->GetElVertices (ElementId (BND, elnr));
        v.Col (0).Range (0, D) = ma->GetPoint<D> (tent->vertex);
        v (D, 0) = tent->tbot;
        for (int n = 0; n < D; n++)
          {
            v.Col (n + 1).Range (0, D) = ma->GetPoint<D> (sel_verts[n]);
            v (D, n + 1) = tent->vertex == sel_verts[n]
                               ? tent->ttop
                               : tent->nbtime[tent->nbv.Pos (sel_verts[n])];
          }
      }
    else // top or bot of tent
      {
        INT<D + 1> vnr = ma->GetElVertices (elnr);
        // determine linear basis function coeffs to use for tent face
        for (int ivert = 0; ivert < vnr.Size (); ivert++)
          {
            v.Col (ivert).Range (0, D) = ma->GetPoint<D> (vnr[ivert]);
            if (vnr[ivert] == tent->vertex)
              v (D, ivert) = top == 1 ? tent->ttop : tent->tbot;
            else
              for (int k = 0; k < tent->nbv.Size (); k++)
                if (vnr[ivert] == tent->nbv[k])
                  v (D, ivert) = tent->nbtime[k];
          }
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
  Matrix<> MakeWavefront (int order, shared_ptr<MeshAccess> ma, double time,
                          shared_ptr<CoefficientFunction> bddatum)
  {
    LocalHeap lh (10000000);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    Matrix<> wavefront (ma->GetNE (), snip * (D + 2));
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D + 1> smir (
            sir, ma->GetTrafo (elnr, lh), -1, lh);
        SIMD_MappedIntegrationRule<D, D> smir_fix (
            sir, ma->GetTrafo (elnr, lh), lh);
        for (int imip = 0; imip < sir.Size (); imip++)
          {
            smir[imip].Point ().Range (0, D)
                = smir_fix[imip].Point ().Range (0, D);
            smir[imip].Point ()[D] = time;
          }
        FlatMatrix<SIMD<double>> bdeval (D + 2, smir.Size (), lh);
        bdeval = 0;
        bddatum->Evaluate (smir, bdeval);
        for (int imip = 0; imip < snip; imip++)
          {
            wavefront (elnr, imip) = bdeval (0, imip / nsimd)[imip % nsimd];
            for (int d = 0; d < D + 1; d++)
              wavefront (elnr, snip + d * snip + imip)
                  = bdeval (d + 1, imip / nsimd)[imip % nsimd];
          }
      }
    return wavefront;
  }

  template <int D>
  double Error (int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront,
                Matrix<> wavefront_corr)
  {
    LocalHeap lh (100000000);
    double error = 0;
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);
        for (int imip = 0; imip < snip; imip++)
          {
            // l2error +=
            // (wavefront(elnr,imip)-wavefront_corr(elnr,imip))*(wavefront(elnr,imip)-wavefront_corr(elnr,imip))*ir[imip].Weight();
            for (int d = 0; d < D + 1; d++)
              error
                  += pow (wavefront (elnr, snip + d * snip + imip)
                              - wavefront_corr (elnr, snip + d * snip + imip),
                          2)
                     * smir[imip / nsimd].GetWeight ()[imip % nsimd];
          }
      }
    return sqrt (error);
  }

  template <int D>
  double Energy (int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront,
                 double wavenumber)
  {
    double energy = 0;
    LocalHeap lh (10000000);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    IntegrationRule ir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = ir.Size ()
               + (ir.Size () % nsimd == 0 ? 0 : nsimd - ir.Size () % nsimd);
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      for (int imip = 0; imip < ir.Size (); imip++)
        for (int d = 0; d < D + 1; d++)
          energy += 0.5
                    * (((d == D ? 1.0 : 0.0) / pow (wavenumber, 2))
                       * wavefront (elnr, snip + d * snip + imip)
                       * wavefront (elnr, snip + d * snip + imip))
                    * ir[imip].Weight ();

    return energy;
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

  template <int D>
  double TentAdiam (Tent *tent, shared_ptr<MeshAccess> ma, double wavespeed)
  {
    double anisotropicdiam = 0;
    int vnumber = tent->nbv.Size ();

    // Array<int> verts(vnumber);
    // verts.Range(2,vnumber) = tent->nbv;
    // verts[0] = tent->vertex;
    // verts[1] = tent->vertex;

    // Array<int> vtime(vnumber);
    // vtime.Range(2,vnumber) = tent->nbtime;
    // vtime[0] = tent->tbot;
    // vtime[1] = tent->ttop;

    for (int k = 0; k < vnumber; k++)
      {
        Vec<D> v1 = ma->GetPoint<D> (tent->vertex);
        Vec<D> v2 = ma->GetPoint<D> (tent->nbv[k]);
        anisotropicdiam = max (
            anisotropicdiam,
            sqrt (L2Norm2 (v1 - v2)
                  + pow (wavespeed * (tent->ttop - tent->nbtime[k]), 2)));
        anisotropicdiam = max (
            anisotropicdiam,
            sqrt (L2Norm2 (v1 - v2)
                  + pow (wavespeed * (tent->tbot - tent->nbtime[k]), 2)));
        for (int j = 0; j < vnumber; j++)
          {
            v1 = ma->GetPoint<D> (tent->nbv[j]);
            v2 = ma->GetPoint<D> (tent->nbv[k]);
            anisotropicdiam = max (
                anisotropicdiam,
                sqrt (L2Norm2 (v1 - v2)
                      + pow (wavespeed * (tent->nbtime[j] - tent->nbtime[k]),
                             2)));
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
             Matrix<> wavefront, double timeshift,
             shared_ptr<CoefficientFunction> bddatum)
             -> Matrix<> //-> shared_ptr<MeshAccess>
         {
           int D = ma->GetDimension ();
           if (D == 1)
             EvolveTents<1> (order, ma, wavespeed, dt, wavefront, timeshift,
                             bddatum);
           else if (D == 2)
             EvolveTents<2> (order, ma, wavespeed, dt, wavefront, timeshift,
                             bddatum);
           else if (D == 3)
             EvolveTents<3> (order, ma, wavespeed, dt, wavefront, timeshift,
                             bddatum);
           return wavefront;
         } //, py::call_guard<py::gil_scoped_release>()
           // py::arg("oder"),py::arg("ma"),py::arg("ws"),py::arg("finaltime"),
         // py::arg("wavefront"), py::arg("timeshift"), py::arg("solname") = ""
  );

  m.def ("EvolveTents",
         [] (int order, shared_ptr<MeshAccess> ma, Vector<double> wavespeed,
             double dt, Matrix<> wavefront, double timeshift,
             shared_ptr<CoefficientFunction> bddatum) -> Matrix<> {
           int D = ma->GetDimension ();
           if (D == 1)
             EvolveTents<1> (order, ma, wavespeed, dt, wavefront, timeshift,
                             bddatum);
           else if (D == 2)
             EvolveTents<2> (order, ma, wavespeed, dt, wavefront, timeshift,
                             bddatum);
           else if (D == 3)
             EvolveTents<3> (order, ma, wavespeed, dt, wavefront, timeshift,
                             bddatum);
           return wavefront;
         });

  m.def ("EvolveTentsMakeWavefront",
         [] (int order, shared_ptr<MeshAccess> ma, double time,
             shared_ptr<CoefficientFunction> bddatum) -> Matrix<> {
           int D = ma->GetDimension ();
           Matrix<> wavefront;
           if (D == 1)
             wavefront = MakeWavefront<1> (order, ma, time, bddatum);
           else if (D == 2)
             wavefront = MakeWavefront<2> (order, ma, time, bddatum);
           else if (D == 3)
             wavefront = MakeWavefront<3> (order, ma, time, bddatum);
           return wavefront;
         });

  m.def ("EvolveTentsError",
         [] (int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront,
             Matrix<> wavefront_corr) -> double {
           int D = ma->GetDimension ();
           double error;
           if (D == 1)
             error = Error<1> (order, ma, wavefront, wavefront_corr);
           else if (D == 2)
             error = Error<2> (order, ma, wavefront, wavefront_corr);
           else if (D == 3)
             error = Error<3> (order, ma, wavefront, wavefront_corr);
           return error;
         });

  m.def ("EvolveTentsEnergy",
         [] (int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront,
             double wavenumber) -> double {
           int D = ma->GetDimension ();
           double energy;
           if (D == 1)
             energy = Energy<1> (order, ma, wavefront, wavenumber);
           else if (D == 2)
             energy = Energy<2> (order, ma, wavefront, wavenumber);
           else if (D == 3)
             energy = Energy<3> (order, ma, wavefront, wavenumber);
           return energy;
         });

  m.def ("EvolveTentsDofs",
         [] (int order, shared_ptr<MeshAccess> ma, double wavespeed,
             double dt) -> int {
           int D = ma->GetDimension ();
           int dofis;
           if (D == 1)
             {
               TentPitchedSlab<1> tps = TentPitchedSlab<1> (ma);
               tps.PitchTents (dt, wavespeed + 1);
               TrefftzWaveFE<2> tel (order, wavespeed);
               dofis = tps.tents.Size () * tel.GetNBasis ();
             }
           else if (D == 2)
             {
               TentPitchedSlab<2> tps = TentPitchedSlab<2> (ma);
               tps.PitchTents (dt, wavespeed + 1);
               TrefftzWaveFE<3> tel (order, wavespeed);
               dofis = tps.tents.Size () * tel.GetNBasis ();
             }
           else if (D == 3)
             {
               TentPitchedSlab<3> tps = TentPitchedSlab<3> (ma);
               tps.PitchTents (dt, wavespeed + 1);
               TrefftzWaveFE<4> tel (order, wavespeed);
               dofis = tps.tents.Size () * tel.GetNBasis ();
             }
           return dofis;
         });

  m.def (
      "EvolveTentsAdiam",
      [] (shared_ptr<MeshAccess> ma, double wavespeed, double dt) -> double {
        int D = ma->GetDimension ();
        double h = 0.0;
        if (D == 1)
          {
            TentPitchedSlab<1> tps = TentPitchedSlab<1> (ma);
            tps.PitchTents (dt, wavespeed + 1);
            RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
              Tent *tent = tps.tents[tentnr];
              h = max (h, TentAdiam<1> (tent, ma, wavespeed));
            });
          }
        else if (D == 2)
          {
            TentPitchedSlab<2> tps = TentPitchedSlab<2> (ma);
            tps.PitchTents (dt, wavespeed + 1);
            RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
              Tent *tent = tps.tents[tentnr];
              h = max (h, TentAdiam<2> (tent, ma, wavespeed));
            });
          }
        else if (D == 3)
          {
            TentPitchedSlab<3> tps = TentPitchedSlab<3> (ma);
            tps.PitchTents (dt, wavespeed + 1);
            RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
              Tent *tent = tps.tents[tentnr];
              h = max (h, TentAdiam<3> (tent, ma, wavespeed));
            });
          }
        return h;
      });
}
#endif // NGS_PYTHON
