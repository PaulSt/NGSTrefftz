#include "evolvetent.hpp"

namespace ngcomp
{
  template <int D>
  inline void
  WaveTents<D>::LapackSolve (SliceMatrix<double> a, SliceVector<double> b)
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

  template <int D> void WaveTents<D>::EvolveTents (double dt)
  {
    // int nthreads = (task_manager) ? task_manager->GetNumThreads() : 1;
    LocalHeap lh (1000 * 1000 * 100, "trefftz tents", 1);

    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    const int nsimd = SIMD<double>::Size ();
    SIMD_IntegrationRule sir (eltyp, order * 2);
    const int snip = sir.Size () * nsimd;

    const int ndomains = ma->GetNDomains ();
    double max_wavespeed = wavespeed[0];
    for (double c : wavespeed)
      max_wavespeed = max (c, max_wavespeed);

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (ma); // collection of tents in timeslab
    tps.PitchTents (
        dt, this->wavespeedcf + make_shared<ConstantCoefficientFunction> (1),
        lh); // adt = time slab height, wavespeed

    cout << "solving " << tps.tents.Size () << " tents ";
    static Timer ttent ("tent", 2);
    static Timer ttentel ("tentel", 2);
    static Timer ttentbnd ("tentbnd", 2);
    static Timer ttentmacro ("tentmacro", 2);
    static Timer ttenteval ("tenteval", 2);

    TrefftzWaveBasis<D>::getInstance ().CreateTB (order);

    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      LocalHeap slh = lh.Split (); // split to threads
      Tent *tent = tps.tents[tentnr];

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      TrefftzWaveFE<D> tel (
          order, max_wavespeed, center,
          TentAdiam (tent)); // TODO: fix scaling with correct wavespeed
      int nbasis = tel.GetNDof ();

      std::unordered_map<int, int> macroel;
      int ndomains = MakeMacroEl (tent->els, macroel);

      FlatMatrix<> elmat (ndomains * nbasis, slh);
      FlatVector<> elvec (ndomains * nbasis, slh);
      elmat = 0;
      elvec = 0;

      Array<FlatMatrix<SIMD<double>>> topdshapes (tent->els.Size ());
      for (auto &tds : topdshapes)
        tds.AssignMemory ((D + 1) * nbasis, sir.Size (), slh);

      // Integrate top and bottom space-like tent faces
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          tel.SetWavespeed (wavespeed[tent->els[elnr]]);
          int eli = ndomains > 1 ? macroel[tent->els[elnr]] : 0;
          SliceMatrix<> subm = elmat.Cols (eli * nbasis, (eli + 1) * nbasis)
                                   .Rows (eli * nbasis, (eli + 1) * nbasis);
          SliceVector<> subv = elvec.Range (eli * nbasis, (eli + 1) * nbasis);
          CalcTentEl (tent->els[elnr], tent, tel, sir, slh, subm, subv,
                      topdshapes[elnr]);
        }

      for (auto fnr : tent->edges)
        {
          Array<int> elnums;
          ma->GetFacetElements (fnr, elnums);
          Array<int> selnums;
          // ma->GetFacetSurfaceElements (fnr, selnums);
          if (elnums.Size () == 1)
            switch (D)
              {
              case 1:
                selnums = ma->GetVertexSurfaceElements (fnr);
                break;
              case 2:
                ma->GetEdgeSurfaceElements (fnr, selnums);
                break;
              case 3:
                {
                  for (int i : Range (ma->GetNSE ()))
                    {
                      auto sel = ElementId (
                          BND, i); //<-- ist das Oberflächenelement
                      auto fnums = ma->GetElFacets (sel);
                      if (fnr == fnums[0]) //<-- das ist dann die Facet-Nr.,
                                           //also unter allen Facets im Mesh,
                                           //kannst du dir wo speichern
                        {
                          selnums.Append (i);
                          break;
                        }
                    }
                }
              }

          // Integrate boundary tent
          if (elnums.Size () == 1 && selnums.Size () == 1)
            {
              tel.SetWavespeed (wavespeed[elnums[0]]);
              int eli = ndomains > 1 ? macroel[elnums[0]] : 0;

              SliceMatrix<> subm
                  = elmat.Cols (eli * nbasis, (eli + 1) * nbasis)
                        .Rows (eli * nbasis, (eli + 1) * nbasis);
              SliceVector<> subv
                  = elvec.Range (eli * nbasis, (eli + 1) * nbasis);
              CalcTentBndEl (selnums[0], tent, tel, sir, slh, subm, subv);
            }

          // Integrate macro bnd inside tent
          else if (elnums.Size () == 2 && ndomains > 1
                   && macroel[elnums[0]] != macroel[elnums[1]])
            {
              CalcTentMacroEl (fnr, elnums, macroel, tent, tel, sir, slh,
                               elmat, elvec);
            }
        }

      // solve
      LapackSolve (elmat, elvec);
      FlatVector<> sol (ndomains * nbasis, &elvec (0));

      // eval solution on top of tent
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          tel.SetWavespeed (wavespeed[tent->els[elnr]]);
          int eli = ndomains > 1 ? macroel[tent->els[elnr]] : 0;
          CalcTentElEval (tent->els[elnr], tent, tel, sir, slh,
                          sol.Range (eli * nbasis, (eli + 1) * nbasis),
                          topdshapes[elnr]);
        }
    }); // end loop over tents
    cout << "solved from " << timeshift;
    timeshift += dt;
    cout << " to " << timeshift << endl;
    // return wavefront;
  }

  template <int D>
  void WaveTents<D>::CalcTentEl (int elnr, Tent *tent,
                                 ScalarMappedElement<D + 1> &tel,
                                 SIMD_IntegrationRule &sir, LocalHeap &slh,
                                 SliceMatrix<> elmat, SliceVector<> elvec,
                                 SliceMatrix<SIMD<double>> simddshapes)
  {
    static Timer tint1 ("tent int calcshape", 2);
    static Timer tint2 ("tent int mat&vec", 2);
    static Timer tint3 ("tent mat*vec", 2);

    HeapReset hr (slh);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    double wavespeed = tel.GetWavespeed ();
    int nbasis = tel.GetNDof ();
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
    vert = TentFaceVerts (tent, elnr, -1);
    linbasis = vert.Row (D);
    faceint.Evaluate (sir, linbasis, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
    tel.CalcShape (smir, simdshapes);
    tel.CalcDShape (smir, simddshapes);
    FlatMatrix<> bbmat (nbasis, (D + 1) * snip, &simddshapes (0, 0)[0]);

    TentDmat (Dmat, vert, -1, wavespeed);

    Vec<D + 1> nnn = TentFaceNormal (vert, 1);

    static bool displaytenterror = true;
    if (L2Norm (nnn.Range (0, D)) * (wavespeed) > nnn[D] && displaytenterror)
      {
        cout << endl << "some tents pitched too high" << endl;
        displaytenterror = false;
      }
    // if(TentFaceArea<D>(vert)<smir_fix[0].GetMeasure()[0]&&
    // tent->nbtime[0]==0)
    //{
    // cout << "tent area problem "<<endl<< TentFaceArea<D>(vert)<< " smaller
    // thamn " <<smir_fix[0].GetMeasure()[0]<< " ratio " <<
    // smir_fix[0].GetJacobiDet()[0]/TentFaceArea<D>(vert) << " vert " << endl;
    // for(int iii=0;iii<D+1;iii++) cout<<vert.Row(iii)<<endl;
    // }

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
          *= sqrt (TentFaceArea (vert) * sir[imip].Weight ());
    AddABt (simdshapes, simdshapes, elmat);
    for (int imip = 0; imip < sir.Size (); imip++)
      simdshapes.Col (imip)
          *= sqrt (TentFaceArea (vert) * sir[imip].Weight ());
    FlatMatrix<> shapes (nbasis, snip, &simdshapes (0, 0)[0]);
    elvec += shapes * wavefront.Row (elnr).Range (0, snip);

    /// Integration over top of tent
    vert = TentFaceVerts (tent, elnr, 1);
    linbasis = vert.Row (D);
    faceint.Evaluate (sir, linbasis, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    // FlatMatrix<SIMD<double>> simddshapes((D+1)*nbasis,sir.Size(),slh);
    tint1.Start ();
    tel.CalcDShape (smir, simddshapes);
    tint1.Stop ();

    TentDmat (Dmat, vert, 1, wavespeed);
    FlatMatrix<> bdbmat ((D + 1) * snip, nbasis, slh);
    bdbmat = 0;
    tint2.Start ();
    for (int imip = 0; imip < sir.Size (); imip++)
      for (int si = 0; si < nsimd; si++)
        for (int r = 0; r < (D + 1); r++)
          for (int d = 0; d < D + 1; d++)
            bdbmat.Row (r * snip + imip * nsimd + si)
                += Dmat (r, d) * sir[imip].Weight ()[si]
                   * bbmat.Col (d * snip + imip * nsimd + si);
    tint2.Stop ();
    tint3.Start ();
    elmat += bbmat * bdbmat;
    tint3.Stop ();
    tint3.AddFlops (bbmat.Height () * bbmat.Width () * bdbmat.Width ());
    // cout << "height " << bbmat.Height() << " width " << bbmat.Width()  <<
    // endl;
  }

  template <int D>
  void WaveTents<D>::CalcTentBndEl (int surfel, Tent *tent,
                                    ScalarMappedElement<D + 1> &tel,
                                    SIMD_IntegrationRule &sir, LocalHeap &slh,
                                    SliceMatrix<> elmat, SliceVector<> elvec)
  {
    HeapReset hr (slh);
    double wavespeed = tel.GetWavespeed ();
    int nbasis = tel.GetNDof ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;

    // get vertices of tent face
    Mat<D + 1> vert = TentFaceVerts (tent, surfel, 0);
    // build normal vector
    Vec<D + 1> n;
    n = -TentFaceNormal (vert, 0);

    if (D == 1) // D=1 special case
      n[0] = sgn_nozero<int> (tent->vertex - tent->nbv[0]);
    n[D] = 0;
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
    Dmat.Row (D).Range (0, D) = -TentFaceArea (vert) * n.Range (0, D);
    FlatMatrix<double> bdbmat ((D + 1) * snip, nbasis, slh);
    bdbmat = 0;

    for (int imip = 0; imip < smir.Size (); imip++)
      smir[imip].Point ()[D] += timeshift;
    FlatMatrix<SIMD<double>> bdeval ((D + 2), sir.Size (), slh);
    bdeval = 0;
    bddatum->Evaluate (smir, bdeval);

    FlatVector<> bdbvec ((D + 1) * snip, slh);
    bdbvec = 0;
    if (ma->GetMaterial (ElementId (BND, surfel)) == "neumann")
      {
        double beta = 0.5;
        double A = TentFaceArea (vert);
        for (int imip = 0; imip < snip; imip++)
          for (int r = 0; r < D; r++)
            for (int d = 0; d < (D + 1); d++)
              bdbmat.Row (r * snip + imip)
                  += (d < D ? -n (d) * beta : 1.0) * (-n (r))
                     * sir[imip / nsimd].Weight ()[imip % nsimd] * A
                     * bbmat.Col (d * snip + imip);
        elmat += bbmat * bdbmat;

        for (int imip = 0; imip < snip; imip++)
          for (int r = 0; r < D; r++)
            for (int d = 0; d < D + 1; d++)
              bdbvec (d * snip + imip)
                  += (d < D ? -n (d) * beta : -1.0) * (-n (r))
                     * bdeval (r + 1, imip / nsimd)[imip % nsimd]
                     * sir[imip / nsimd].Weight ()[imip % nsimd] * A;
        elvec += bbmat * bdbvec;
      }
    else
      { // dirichlet
        // double alpha = 0.5;
        // Dmat(D,D) = TentFaceArea(vert)*alpha;

        FlatMatrix<SIMD<double>> wavespeed (1, sir.Size (), slh);
        auto localwavespeedcf = make_shared<ConstantCoefficientFunction> (1)
                                / (this->wavespeedcf);
        localwavespeedcf->Evaluate (smir, wavespeed);
        double DmatDD = TentFaceArea (vert);

        for (int imip = 0; imip < snip; imip++)
          // for(int r=0;r<(D+1);r++) // r=D since only last row non-zero
          for (int d = 0; d < D + 1; d++)
            {
              Dmat (D, D) = DmatDD * wavespeed (0, imip / nsimd)[imip % nsimd];
              bdbmat.Row (D * snip + imip)
                  += Dmat (D, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                     * bbmat.Col (d * snip + imip);
            }
        elmat += bbmat * bdbmat;

        Dmat (D, D) *= -1.0;
        for (int imip = 0; imip < snip; imip++)
          for (int r = 0; r < (D + 1); r++)
            {
              Dmat (D, D)
                  = -DmatDD * wavespeed (0, imip / nsimd)[imip % nsimd];
              bdbvec (r * snip + imip)
                  += Dmat (D, r) * sir[imip / nsimd].Weight ()[imip % nsimd]
                     * bdeval (
                         D + 1,
                         imip / nsimd)[imip % nsimd]; // use Dmat transposed
            }
        elvec -= bbmat * bdbvec;
      }
  }

  template <int D>
  void
  WaveTents<D>::CalcTentMacroEl (int fnr, const Array<int> &elnums,
                                 std::unordered_map<int, int> &macroel,
                                 Tent *tent, TrefftzWaveFE<D> &tel,
                                 SIMD_IntegrationRule &sir, LocalHeap &slh,
                                 SliceMatrix<> elmat, SliceVector<> elvec)
  {
    int nbasis = tel.GetNDof ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;

    Array<int> fnums;
    Array<int> orient;
    switch (D)
      {
      case 1:
        fnums.Append (fnr);
        orient.Append (1);
        break;
      case 2:
        ma->GetElEdges (elnums[0], fnums, orient);
        break;
      case 3:
        ma->GetElFaces (elnums[0], fnums, orient);
        break;
      }

    // get vertices of tent face
    // Mat<D+1> vert = TentFaceVerts(tent, surfel, 0);
    // auto sel_verts = ma->GetElVertices(ElementId(BND,elnr));
    Mat<D + 1, D + 1> vert;
    Array<int> sel_verts (D);
    ma->GetFacetPNums (fnr, sel_verts);
    vert.Col (0).Range (0, D) = ma->GetPoint<D> (tent->vertex);
    vert (D, 0) = tent->tbot;
    for (int n = 0; n < D; n++)
      {
        vert.Col (n + 1).Range (0, D) = ma->GetPoint<D> (sel_verts[n]);
        vert (D, n + 1) = tent->vertex == sel_verts[n]
                              ? tent->ttop
                              : tent->nbtime[tent->nbv.Pos (sel_verts[n])];
      }

    // build normal vector
    Vec<D + 1> n;
    n = -orient[fnums.Pos (fnr)] * TentFaceNormal (vert, 0);
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
    for (int imip = 0; imip < snip; imip++)
      smir[imip].Point ()
          = map * sir[imip].operator Vec<D, SIMD<double>> () + shift;

    FlatMatrix<> *bbmat[2];

    tel.SetWavespeed (this->wavespeed[elnums[0]]);
    FlatMatrix<SIMD<double>> simddshapes1 ((D + 1) * nbasis, sir.Size (), slh);
    tel.CalcDShape (smir, simddshapes1);
    bbmat[0]
        = new FlatMatrix<> (nbasis, (D + 1) * snip, &simddshapes1 (0, 0)[0]);

    tel.SetWavespeed (this->wavespeed[elnums[1]]);
    FlatMatrix<SIMD<double>> simddshapes2 ((D + 1) * nbasis, sir.Size (), slh);
    tel.CalcDShape (smir, simddshapes2);
    bbmat[1]
        = new FlatMatrix<> (nbasis, (D + 1) * snip, &simddshapes2 (0, 0)[0]);

    Mat<D + 1> Dmat1 = 0;
    Dmat1.Row (D).Range (0, D)
        = -0.5 * TentFaceArea (vert) * n.Range (0, D); // 0.5 for DG average
    Dmat1.Col (D).Range (0, D) = -0.5 * TentFaceArea (vert) * n.Range (0, D);

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
          for (int el = 0; el < 4; el++)
            {
              bdbmat[el]->Row (r * snip + imip)
                  += pow (-1, el / 2) * Dmat1 (r, d)
                     * sir[imip / nsimd].Weight ()[imip % nsimd]
                     * bbmat[el % 2]->Col (d * snip + imip);
            }

    for (int el = 0; el < 4; el++)
      {
        int in = macroel[elnums[el / 2]];
        int out = macroel[elnums[el % 2]];
        elmat.Cols (out * nbasis, (out + 1) * nbasis)
            .Rows (in * nbasis, (in + 1) * nbasis)
            += *bbmat[el / 2] * (*bdbmat[el]);
      }
  }

  template <int D>
  void WaveTents<D>::CalcTentElEval (int elnr, Tent *tent,
                                     ScalarMappedElement<D + 1> &tel,
                                     SIMD_IntegrationRule &sir, LocalHeap &slh,
                                     SliceVector<> sol,
                                     SliceMatrix<SIMD<double>> simddshapes)
  {
    HeapReset hr (slh);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    int nbasis = tel.GetNDof ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces

    SIMD_MappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (elnr, slh),
                                               -1, slh);
    SIMD_MappedIntegrationRule<D, D> smir_fix (sir, ma->GetTrafo (elnr, slh),
                                               slh);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ().Range (0, D) = smir_fix[imip].Point ().Range (0, D);

    Mat<D + 1> v = TentFaceVerts (tent, elnr, 1);
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
  void WaveTents<D>::TentDmat (Mat<D + 1> &Dmat, Mat<D + 1> v, int top,
                               double wavespeed)
  {
    Vec<D + 1> n = TentFaceNormal (v, top);
    Dmat = n (D) * Id<D + 1> ();
    Dmat.Row (D).Range (0, D) = -n.Range (0, D);
    Dmat.Col (D).Range (0, D) = -n.Range (0, D);
    Dmat (D, D) *= 1.0 / (wavespeed * wavespeed);
    Dmat *= TentFaceArea (v);
  }

  // returns matrix where cols correspond to vertex coordinates of the
  // space-time element
  template <int D>
  Mat<D + 1, D + 1> WaveTents<D>::TentFaceVerts (Tent *tent, int elnr, int top)
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
        INT<D + 1> vnr = ma->GetElVertices (ElementId (VOL, elnr));
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

  template <int D> double WaveTents<D>::TentFaceArea (Mat<D + 1, D + 1> ve)
  {
    switch (D)
      {
      case 1:
        return L2Norm (ve.Col (0) - ve.Col (1));
        break;
      case 2:
        {
          // double a = L2Norm(ve.Col(0)-ve.Col(1));
          // double b = L2Norm(ve.Col(1)-ve.Col(2));
          // double c = L2Norm(ve.Col(0)-ve.Col(2));
          // SwapIfGreater<>(a,b); SwapIfGreater<>(a,c); SwapIfGreater<>(b,c);
          // return 0.25 * sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)))  *2.0;
          Vec<3> ba = ve.Col (0) - ve.Col (1);
          Vec<3> baa = ve.Col (0) - ve.Col (2);
          return L2Norm (Cross (ba, baa)) / 2.0 * 2.0;
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
                 / (192.0 * u * v * w) * 6.0;
        }
      }
  }

  template <int D>
  Vec<D + 1> WaveTents<D>::TentFaceNormal (Mat<D + 1, D + 1> v, int top)
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
          for (int d = 1; d < D + 1; d++)
            v.Col (d) = v.Col (0) - v.Col (d);

          for (unsigned int i = 0; i < D + 1; i++)
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
    if (top == 1)
      normv *= sgn_nozero<double> (normv[D]);
    else if (top == -1)
      normv *= (-sgn_nozero<double> (normv[D]));
    else if (top == 0)
      normv[D] = 0; // time-like faces only
    normv /= L2Norm (normv);
    return normv;
  }

  template <int D>
  Matrix<>
  WaveTents<D>::MakeWavefront (shared_ptr<CoefficientFunction> bddatum,
                               double time)
  {
    LocalHeap lh (1000 * 1000 * 100, "make wavefront", 1);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    Matrix<> wf (ma->GetNE (), snip * (D + 2));
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
            wf (elnr, imip) = bdeval (0, imip / nsimd)[imip % nsimd];
            for (int d = 0; d < D + 1; d++)
              wf (elnr, snip + d * snip + imip)
                  = bdeval (d + 1, imip / nsimd)[imip % nsimd];
          }
      }
    return wf;
  }

  template <int D>
  double WaveTents<D>::Error (Matrix<> wavefront, Matrix<> wavefront_corr)
  {
    LocalHeap lh (1000 * 1000 * 100, "error", 1);
    double error = 0;
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);

        for (int imip = 0; imip < snip; imip++)
          {
            IntegrationRule ir (eltyp, 0);
            MappedIntegrationPoint<D, D> mip (
                ir[0], ma->GetTrafo (ElementId (0), lh));
            for (int d = 0; d < D; d++)
              mip.Point ()[d] = smir[imip / nsimd].Point ()[d][0];
            // mip.Point()[D] = timeshift;
            for (int d = 0; d < D + 1; d++)
              error += pow (wavespeedcf->Evaluate (mip), -2 * (d == D))
                       * pow (
                           wavefront (elnr, snip + d * snip + imip)
                               - wavefront_corr (elnr, snip + d * snip + imip),
                           2)
                       * smir[imip / nsimd].GetWeight ()[imip % nsimd];
            // error +=
            // (pow(wavefront(elnr,snip+d*snip+imip)-wavefront_corr(elnr,snip+d*snip+imip),2)*
            // smir[imip/nsimd].GetWeight()[imip%nsimd]);
          }
      }
    return sqrt (error);
  }

  template <int D>
  double WaveTents<D>::L2Error (Matrix<> wavefront, Matrix<> wavefront_corr)
  {
    LocalHeap lh (1000 * 1000 * 100, "l2error", 1);
    double l2error = 0;
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);
        for (int imip = 0; imip < snip; imip++)
          {
            l2error += (wavefront (elnr, imip) - wavefront_corr (elnr, imip))
                       * (wavefront (elnr, imip) - wavefront_corr (elnr, imip))
                       * smir[imip / nsimd].GetWeight ()[imip % nsimd];
          }
      }
    return sqrt (l2error);
  }

  template <int D> double WaveTents<D>::Energy (Matrix<> wavefront)
  {
    double energy = 0;
    LocalHeap lh (1000 * 1000 * 100, "energy", 1);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    for (int elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);
        FlatMatrix<SIMD<double>> wavespeed (1, smir.Size (), lh);
        wavespeedcf->Evaluate (smir, wavespeed);
        for (int imip = 0; imip < snip; imip++)
          for (int d = 0; d < D + 1; d++)
            energy += 0.5
                      * (pow (wavespeed (0, imip / nsimd)[imip % nsimd],
                              -2 * (d == D))
                         * pow (wavefront (elnr, snip + d * snip + imip), 2)
                         * smir[imip / nsimd].GetWeight ()[imip % nsimd]);
      }

    return energy;
  }

  template <int D>
  template <typename T>
  void WaveTents<D>::SwapIfGreater (T &a, T &b)
  {
    if (a < b)
      {
        T tmp (a);
        a = b;
        b = tmp;
      }
  }

  template <int D> double WaveTents<D>::TentAdiam (Tent *tent)
  {
    LocalHeap lh (1000 * 1000 * 100);
    int vnumber = tent->nbv.Size ();

    // Array<int> verts(vnumber);
    // verts.Range(2,vnumber) = tent->nbv;
    // verts[0] = tent->vertex;
    // verts[1] = tent->vertex;

    // Array<int> vtime(vnumber);
    // vtime.Range(2,vnumber) = tent->nbtime;
    // vtime[0] = tent->tbot;
    // vtime[1] = tent->ttop;

    double c = wavespeed[tent->els[0]];
    // for(auto el : tent->els) c = max(c,wavespeed[el]);
    IntegrationRule ir (ma->GetElType (ElementId (VOL, 0)), 0);
    ElementTransformation &trafo = ma->GetTrafo (0, lh);
    MappedIntegrationPoint<D, D> mip (ir[0], trafo);
    Vec<D> v1 = ma->GetPoint<D> (tent->vertex);
    double c1 = wavespeedcf->Evaluate (mip);
    double anisotropicdiam = c1 * tent->ttop - c1 * tent->tbot;

    for (int k = 0; k < vnumber; k++)
      {
        Vec<D> v1 = ma->GetPoint<D> (tent->vertex);
        Vec<D> v2 = ma->GetPoint<D> (tent->nbv[k]);
        mip.Point () = v1;
        double c1 = wavespeedcf->Evaluate (mip);
        mip.Point () = v2;
        double c2 = wavespeedcf->Evaluate (mip);

        anisotropicdiam
            = max (anisotropicdiam,
                   sqrt (L2Norm2 (v1 - v2)
                         + pow (c1 * tent->ttop - c2 * tent->nbtime[k], 2)));
        anisotropicdiam
            = max (anisotropicdiam,
                   sqrt (L2Norm2 (v1 - v2)
                         + pow (c1 * tent->tbot - c2 * tent->nbtime[k], 2)));
        for (int j = 0; j < vnumber; j++)
          {
            v1 = ma->GetPoint<D> (tent->nbv[j]);
            mip.Point () = v1;
            c1 = wavespeedcf->Evaluate (mip);
            anisotropicdiam = max (
                anisotropicdiam,
                sqrt (L2Norm2 (v1 - v2)
                      + pow (c1 * tent->nbtime[j] - c2 * tent->nbtime[k], 2)));
          }
      }

    return anisotropicdiam;
  }

  template <int D> double WaveTents<D>::MaxAdiam (double dt)
  {
    double h = 0.0;
    TentPitchedSlab<D> tps = TentPitchedSlab<D> (ma);
    LocalHeap lh (1000 * 1000 * 100);
    tps.PitchTents (
        dt, this->wavespeedcf + make_shared<ConstantCoefficientFunction> (1),
        lh);
    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      Tent *tent = tps.tents[tentnr];
      h = max (h, TentAdiam (tent));
    });
    return h;
  }

  template <int D>
  inline int WaveTents<D>::MakeMacroEl (const Array<int> &tentel,
                                        std::unordered_map<int, int> &macroel)
  {
    // TODO fix if macro elements do not share faces
    int nrmacroel = 0;
    for (int i = 0; i < tentel.Size (); i++)
      {
        int j = 0;
        while (wavespeed[tentel[i]] != wavespeed[tentel[j]])
          j++;
        if (j == i)
          macroel[tentel[i]] = nrmacroel++;
        else
          macroel[tentel[i]] = macroel[tentel[j]];
      }
    return nrmacroel;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D> void GppwTents<D>::EvolveTents (double dt)
  {
    LocalHeap lh (1000 * 1000 * 100, "gppw tents", 1);

    shared_ptr<MeshAccess> ma = this->ma;
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    const int nsimd = SIMD<double>::Size ();
    SIMD_IntegrationRule sir (eltyp, this->order * 2);
    const int snip = sir.Size () * nsimd;

    TentPitchedSlab<D> tps
        = TentPitchedSlab<D> (this->ma); // collection of tents in timeslab
    tps.PitchTents (
        dt, this->wavespeedcf + make_shared<ConstantCoefficientFunction> (1),
        lh);

    cout << "solving gppw " << tps.tents.Size () << " tents " << endl;

    RunParallelDependency (tps.tent_dependency, [&] (int tentnr) {
      LocalHeap slh = lh.Split (); // split to threads
      Tent *tent = tps.tents[tentnr];

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      double tentsize = TentAdiam (tent);

      TrefftzGppwFE<D> tel (this->gamma[tent->vertex], this->order, center,
                            tentsize);
      int nbasis = tel.GetNDof ();

      FlatMatrix<> elmat (nbasis, slh);
      FlatVector<> elvec (nbasis, slh);
      elmat = 0;
      elvec = 0;

      Array<FlatMatrix<SIMD<double>>> topdshapes (tent->els.Size ());
      for (auto &tds : topdshapes)
        tds.AssignMemory ((D + 1) * nbasis, sir.Size (), slh);

      // Integrate top and bottom space-like tent faces
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          CalcTentEl (tent->els[elnr], tent, tel, sir, slh, elmat, elvec,
                      topdshapes[elnr]);
        }

      for (auto fnr : tent->edges)
        {
          Array<int> elnums;
          ma->GetFacetElements (fnr, elnums);
          Array<int> selnums;
          // ma->GetFacetSurfaceElements (fnr, selnums);
          if (elnums.Size () == 1)
            switch (D)
              {
              case 1:
                selnums = ma->GetVertexSurfaceElements (fnr);
                break;
              case 2:
                ma->GetEdgeSurfaceElements (fnr, selnums);
                break;
              case 3:
                cout << "not impl" << endl;
                break;
              }

          // Integrate boundary tent
          if (elnums.Size () == 1 && selnums.Size () == 1)
            this->CalcTentBndEl (selnums[0], tent, tel, sir, slh, elmat,
                                 elvec);
        }

      // integrate volume of tent here
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          /// Integration over bot and top volume of tent element
          for (int part = -1; part <= 1; part += 2)
            {
              // if(tent->nbtime[elnr]==(part<0?tent->tbot:tent->ttop))
              // continue;
              HeapReset hr (slh);
              const ELEMENT_TYPE eltyp = (D == 2) ? ET_TET : ET_TRIG;
              SIMD_IntegrationRule vsir (eltyp, this->order * 2);
              int nbasis = tel.GetNDof ();

              Vec<D + 1> shift;
              shift.Range (0, D) = ma->GetPoint<D> (tent->vertex);
              shift[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
              // shift[D] = tent->nbtime[elnr];
              // shift[D] = part==1 ? tent->nbtime[elnr] : tent->tbot;
              Mat<D + 1, D + 1> map
                  = this->TentFaceVerts (tent, tent->els[elnr], part);

              for (int i = 0; i < D + 1; i++)
                map.Col (i) -= shift;
              double vol = abs (
                  Det (map)); // no need for this * (D==2?6.0:2.0); bc of Det
              if (vol < 10e-16)
                continue;

              SIMD_MappedIntegrationRule<D + 1, D + 1> vsmir (
                  vsir, ma->GetTrafo (tent->els[elnr], slh), -1, slh);
              for (int imip = 0; imip < vsir.Size (); imip++)
                vsmir[imip].Point ()
                    = map * vsir[imip].operator Vec<D + 1, SIMD<double>> ()
                      + shift;

              FlatMatrix<SIMD<double>> wavespeed (1, vsir.Size (), slh);
              auto localwavespeedcf
                  = make_shared<ConstantCoefficientFunction> (1)
                    / (this->wavespeedcf * this->wavespeedcf);
              // auto localwavespeedcf = this->wavespeedcf;
              localwavespeedcf->Evaluate (vsmir, wavespeed);

              FlatMatrix<SIMD<double>> simdddshapes ((D + 1) * nbasis,
                                                     vsir.Size (), slh);
              tel.CalcDDSpecialShape (vsmir, simdddshapes, wavespeed);
              for (int imip = 0; imip < vsir.Size (); imip++)
                simdddshapes.Col (imip) *= vol * vsir[imip].Weight ();
              FlatMatrix<SIMD<double>> simdddshapes2 (
                  nbasis, (D + 1) * vsir.Size (), &simdddshapes (0, 0));

              FlatMatrix<SIMD<double>> simddshapes ((D + 1) * nbasis,
                                                    vsir.Size (), slh);
              tel.CalcDShape (vsmir, simddshapes);
              FlatMatrix<SIMD<double>> simddshapes2 (
                  nbasis, (D + 1) * vsir.Size (), &simddshapes (0, 0));

              AddABt (simdddshapes2, simddshapes2, elmat);

              // volume correction term
              double cmax = abs (wavespeed (0, 0)[0]);
              for (int j = 0; j < vsir.Size ();
                   j++) // auto ws : wavespeed.AsVector())
                for (int i = 0; i < nsimd; i++)
                  cmax = max (cmax, abs (wavespeed (0, j)[i]));
              FlatMatrix<SIMD<double>> mu (1, vsir.Size (), slh);
              localwavespeedcf
                  = make_shared<ConstantCoefficientFunction> (tentsize / cmax)
                    * this->wavespeedcf * this->wavespeedcf;
              localwavespeedcf->Evaluate (vsmir, mu);
              FlatMatrix<SIMD<double>> simdddshapescor ((D + 1) * nbasis,
                                                        vsir.Size (), slh);
              tel.CalcDDSpecialShape (vsmir, simdddshapescor, wavespeed, mu);
              AddABt (FlatMatrix<SIMD<double>> (nbasis, (D + 1) * vsir.Size (),
                                                &simdddshapescor (0, 0)),
                      simdddshapes2, elmat);
            }
        }

      // solve
      LapackSolve (elmat, elvec);
      FlatVector<> sol (nbasis, &elvec (0));

      // eval solution on top of tent
      for (int elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          this->CalcTentElEval (tent->els[elnr], tent, tel, sir, slh, sol,
                                topdshapes[elnr]);
        }
    }); // end loop over tents
    cout << "solved from " << this->timeshift;
    this->timeshift += dt;
    cout << " to " << this->timeshift << endl;
  }

  template <int D>
  void GppwTents<D>::CalcTentEl (int elnr, Tent *tent,
                                 ScalarMappedElement<D + 1> &tel,
                                 SIMD_IntegrationRule &sir, LocalHeap &slh,
                                 SliceMatrix<> elmat, SliceVector<> elvec,
                                 SliceMatrix<SIMD<double>> simddshapes)
  {
    HeapReset hr (slh);
    const ELEMENT_TYPE eltyp
        = (D == 3) ? ET_TET : ((D == 2) ? ET_TRIG : ET_SEGM);
    int nbasis = tel.GetNDof ();
    int nsimd = SIMD<double>::Size ();
    int snip = sir.Size () * nsimd;
    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces

    shared_ptr<MeshAccess> ma = this->ma;
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
    vert = this->TentFaceVerts (tent, elnr, -1);
    linbasis = vert.Row (D);
    faceint.Evaluate (sir, linbasis, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
    tel.CalcShape (smir, simdshapes);
    tel.CalcDShape (smir, simddshapes);
    FlatMatrix<> bbmat (nbasis, (D + 1) * snip, &simddshapes (0, 0)[0]);

    this->TentDmat (Dmat, vert, -1, 1);
    double DmatDD = Dmat (D, D);

    FlatMatrix<SIMD<double>> wavespeed (1, sir.Size (), slh);
    auto localwavespeedcf = make_shared<ConstantCoefficientFunction> (1)
                            / (this->wavespeedcf * this->wavespeedcf);
    // auto localwavespeedcf = this->wavespeedcf;
    localwavespeedcf->Evaluate (smir_fix, wavespeed);
    // for(int s=0;s<smir.Size();s++)
    // cout << "point " << smir[s].Point() << " speed " << wavespeed(0,s) <<
    // endl;

    FlatVector<> bdbvec ((D + 1) * snip, slh);
    bdbvec = 0;
    for (int imip = 0; imip < snip; imip++)
      for (int r = 0; r < (D + 1); r++)
        for (int d = 0; d < D + 1; d++)
          {
            Dmat (D, D) = DmatDD * wavespeed (0, imip / nsimd)[imip % nsimd];
            bdbvec (r * snip + imip)
                += Dmat (r, d) * sir[imip / nsimd].Weight ()[imip % nsimd]
                   * this->wavefront (elnr, snip + d * snip + imip);
          }
    elvec -= bbmat * bdbvec;

    // stabilization to recover second order solution
    for (int imip = 0; imip < sir.Size (); imip++)
      simdshapes.Col (imip)
          *= sqrt (this->TentFaceArea (vert) * sir[imip].Weight ());
    AddABt (simdshapes, simdshapes, elmat);
    for (int imip = 0; imip < sir.Size (); imip++)
      simdshapes.Col (imip)
          *= sqrt (this->TentFaceArea (vert) * sir[imip].Weight ());
    FlatMatrix<> shapes (nbasis, snip, &simdshapes (0, 0)[0]);
    elvec += shapes * this->wavefront.Row (elnr).Range (0, snip);

    /// Integration over top of tent
    vert = this->TentFaceVerts (tent, elnr, 1);
    linbasis = vert.Row (D);
    faceint.Evaluate (sir, linbasis, mirtimes);
    for (int imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    tel.CalcDShape (smir, simddshapes);

    this->TentDmat (Dmat, vert, 1, 1);
    DmatDD = Dmat (D, D);
    FlatMatrix<> bdbmat ((D + 1) * snip, nbasis, slh);
    bdbmat = 0;
    for (int imip = 0; imip < sir.Size (); imip++)
      for (int si = 0; si < nsimd; si++)
        for (int r = 0; r < (D + 1); r++)
          for (int d = 0; d < D + 1; d++)
            {
              Dmat (D, D) = DmatDD * wavespeed (0, imip)[si];
              bdbmat.Row (r * snip + imip * nsimd + si)
                  += Dmat (r, d) * sir[imip].Weight ()[si]
                     * bbmat.Col (d * snip + imip * nsimd + si);
            }
    elmat += bbmat * bdbmat;
  }

}

template class WaveTents<1>;
template class WaveTents<2>;
template class WaveTents<3>;

template class GppwTents<1>;
template class GppwTents<2>;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>

template <class T, int D>
void DeclareETClass (py::module &m, std::string typestr)
{
  using PyETclass = T;
  std::string pyclass_name = typestr;
  // py::class_<PyETclass,shared_ptr<PyETclass> >(m, "EvolveTent")
  py::class_<PyETclass, shared_ptr<PyETclass>, TrefftzTents> (
      m, pyclass_name.c_str ()) //, py::buffer_protocol(), py::dynamic_attr())
                                //.def(py::init<>())
      .def ("EvolveTents", &PyETclass::EvolveTents)
      .def ("MakeWavefront", &PyETclass::MakeWavefront)
      .def ("SetWavefront", &PyETclass::SetWavefront)
      .def ("GetWavefront", &PyETclass::GetWavefront)
      .def ("Error", &PyETclass::Error)
      .def ("L2Error", &PyETclass::L2Error)
      .def ("Energy", &PyETclass::Energy)
      .def ("MaxAdiam", &PyETclass::MaxAdiam)
      .def ("LocalDofs", &PyETclass::LocalDofs)
      .def ("NrTents", &PyETclass::NrTents)
      .def ("GetOrder", &PyETclass::GetOrder)
      .def ("GetSpaceDim", &PyETclass::GetSpaceDim)
      .def ("GetInitmesh", &PyETclass::GetInitmesh);
}

void ExportEvolveTent (py::module m)
{
  py::class_<TrefftzTents, shared_ptr<TrefftzTents>> (
      m, "TrefftzTents"); //, py::buffer_protocol(), py::dynamic_attr())

  DeclareETClass<WaveTents<1>, 1> (m, "WaveTents1");
  DeclareETClass<WaveTents<2>, 2> (m, "WaveTents2");
  DeclareETClass<WaveTents<3>, 3> (m, "WaveTents3");

  DeclareETClass<GppwTents<1>, 1> (m, "GppwTents1");
  DeclareETClass<GppwTents<2>, 2> (m, "GppwTents2");

  // m.def("WaveTents", [](int order, shared_ptr<MeshAccess> ma, double
  // wavespeed, shared_ptr<CoefficientFunction> bddatum) ->
  // shared_ptr<TrefftzTents>
  //{
  ////TrefftzTents* nla = new WaveTents<2>(order, ma, wavespeed, bddatum);
  // shared_ptr<TrefftzTents> tr;
  // int D = ma->GetDimension();
  ////return make_shared<WaveTents<2>>(order,ma,wavespeed,bddatum);
  // if(D==1)
  // tr = make_shared<WaveTents<1>>(order,ma,wavespeed,bddatum);
  // else if(D==2)
  // tr = make_shared<WaveTents<2>>(order,ma,wavespeed,bddatum);
  // else if(D==3)
  // tr = make_shared<WaveTents<3>>(order,ma,wavespeed,bddatum);
  // return tr;
  ////return shared_ptr<TrefftzTents>(new WaveTents<2>(order, ma, wavespeed,
  ///bddatum));

  //});

  m.def (
      "WaveTents",
      [] (int order, shared_ptr<MeshAccess> ma,
          shared_ptr<CoefficientFunction> wavespeedcf,
          shared_ptr<CoefficientFunction> bddatum,
          bool QT) -> shared_ptr<TrefftzTents> {
        // TrefftzTents* nla = new WaveTents<2>(order, ma, wavespeed, bddatum);
        shared_ptr<TrefftzTents> tr;
        int D = ma->GetDimension ();
        // return make_shared<WaveTents<2>>(order,ma,wavespeed,bddatum);
        if (!QT)
          {
            if (D == 1)
              tr = make_shared<WaveTents<1>> (order, ma, wavespeedcf, bddatum);
            else if (D == 2)
              tr = make_shared<WaveTents<2>> (order, ma, wavespeedcf, bddatum);
            else if (D == 3)
              tr = make_shared<WaveTents<3>> (order, ma, wavespeedcf, bddatum);
          }
        else
          {
            if (D == 1)
              tr = make_shared<GppwTents<1>> (order, ma, wavespeedcf, bddatum);
            else if (D == 2)
              tr = make_shared<GppwTents<2>> (order, ma, wavespeedcf, bddatum);
          }
        return tr;
        // return shared_ptr<TrefftzTents>(new WaveTents<2>(order, ma,
        // wavespeed, bddatum));
      },
      "Create Wavetent object", py::arg ("order"), py::arg ("ma"),
      py::arg ("wavespeedcf"), py::arg ("bddatum"), py::arg ("QT") = false);

  m.def ("WaveTents",
         [] (int order, shared_ptr<MeshAccess> ma,
             shared_ptr<CoefficientFunction> wavespeedcf,
             shared_ptr<CoefficientFunction> bddatum,
             vector<shared_ptr<CoefficientFunction>> taylorcf)
             -> shared_ptr<TrefftzTents> {
           // TrefftzTents* nla = new WaveTents<2>(order, ma, wavespeed,
           // bddatum);
           shared_ptr<TrefftzTents> tr;
           int D = ma->GetDimension ();
           // return make_shared<WaveTents<2>>(order,ma,wavespeed,bddatum);
           if (D == 1)
             tr = make_shared<GppwTents<1>> (order, ma, wavespeedcf, bddatum,
                                             taylorcf);
           else if (D == 2)
             tr = make_shared<GppwTents<2>> (order, ma, wavespeedcf, bddatum,
                                             taylorcf);
           return tr;
           // return shared_ptr<TrefftzTents>(new WaveTents<2>(order, ma,
           // wavespeed, bddatum));
         });
}
#endif // NGS_PYTHON
