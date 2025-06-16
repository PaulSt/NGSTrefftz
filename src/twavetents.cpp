#include "twavetents.hpp"
#include <h1lofe.hpp>
#include <paralleldepend.hpp>
#include <fem.hpp>
#include "trefftzfespace.hpp"

namespace ngfem
{
  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_STMappedIntegrationRule<DIM_ELEMENT, DIM_SPACE>::Print (
      ostream &ost) const
  {
    ost << "simd-mir, size = " << mips.Size () << endl;
    for (size_t i = 0; i < mips.Size (); i++)
      mips[i].Print (ost);
  }

  template class SIMD_STMappedIntegrationRule<1, 2>;
  template class SIMD_STMappedIntegrationRule<2, 3>;
  template class SIMD_STMappedIntegrationRule<3, 4>;
}

namespace ngcomp
{
  using ngstents::RunParallelDependency;

  template <int D>
  void SetWavespeed (ScalarMappedElement<D> &tel, double wavespeed)
  {
    Vec<D> scale = tel.GetScale ();
    scale[D - 1] = scale[0] * wavespeed;
    tel.SetScale (scale);
  }

  template <typename T> int sgn_nozero (T val)
  {
    return (T (0) <= val) - (val < T (0));
  }

  template <int D>
  inline void TWaveTents<D>::Solve (FlatMatrix<double> a, FlatVector<double> b)
  {
    CalcInverse (a, INVERSE_LIB::INV_LAPACK);
    Vector<> c = a * b;
    b = c;
  }

  template <int D> void TWaveTents<D>::Propagate ()
  {
    // int nthreads = (task_manager) ? task_manager->GetNumThreads() : 1;
    LocalHeap lh (1000 * 1000 * 1000, "trefftz tents", 1);

    SIMD_IntegrationRule sir (eltyp, order * 2);
    // const int ndomains = ma->GetNDomains();
    double max_wavespeed = wavespeed[0];
    for (double c : wavespeed)
      max_wavespeed = max (c, max_wavespeed);

    // cout << "solving " << tps->GetNTents() << " tents ";
    static Timer ttent ("tent");
    static Timer ttentel ("tentel");
    static Timer ttentbnd ("tentbnd");
    static Timer ttentmacro ("tentmacro");
    static Timer ttenteval ("tenteval");

    CSR basismat = TWaveBasis<D + 1>::Basis (order, 0, fosystem);

    RunParallelDependency (tps->tent_dependency, [&] (int tentnr) {
      LocalHeap slh = lh.Split (); // split to threads
      const Tent *tent = &tps->GetTent (tentnr);

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      ScalarMappedElement<D + 1> tel (nbasis, order, basismat, ET_TET, center,
                                      1.0 / TentAdiam (tent));

      std::unordered_map<int, int> macroel;
      int ndomains = MakeMacroEl (tent->els, macroel);

      FlatMatrix<> elmat (ndomains * nbasis, slh);
      FlatVector<> elvec (ndomains * nbasis, slh);
      elmat = 0;
      elvec = 0;

      for (auto fnr : tent->internal_facets)
        {
          Array<int> elnums;
          ma->GetFacetElements (fnr, elnums);
          Array<int> selnums;
          // ma->GetFacetSurfaceElements (fnr, selnums);
          if (elnums.Size () == 1)
            GetFacetSurfaceElement (ma, fnr, selnums);

          // Integrate boundary tent
          if (elnums.Size () == 1 && selnums.Size () == 1)
            {
              SetWavespeed (tel, wavespeed[elnums[0]]);
              int eli = ndomains > 1 ? macroel[elnums[0]] : 0;

              SliceMatrix<> subm
                  = elmat.Cols (eli * nbasis, (eli + 1) * nbasis)
                        .Rows (eli * nbasis, (eli + 1) * nbasis);
              FlatVector<> subv (nbasis, &elvec (eli * nbasis));
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

      Array<FlatMatrix<SIMD<double>>> topdshapes (tent->els.Size ());
      for (auto &tds : topdshapes)
        tds.AssignMemory ((D + 1) * nbasis, sir.Size (), slh);
      // Integrate top and bottom space-like tent faces
      for (size_t elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          SetWavespeed (tel, wavespeed[tent->els[elnr]]);
          int eli = ndomains > 1 ? macroel[tent->els[elnr]] : 0;
          SliceMatrix<> subm = elmat.Cols (eli * nbasis, (eli + 1) * nbasis)
                                   .Rows (eli * nbasis, (eli + 1) * nbasis);
          FlatVector<> subv (nbasis, &elvec (eli * nbasis));
          double bla = wavespeed[tent->els[elnr]];
          CalcTentEl (
              tent->els[elnr], tent, tel, [&] (int) { return bla; }, sir, slh,
              subm, subv, topdshapes[elnr]);
        }

      // solve
      Solve (elmat, elvec);
      FlatVector<> sol (ndomains * nbasis, &elvec (0));

      // eval solution on top of tent
      for (size_t elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          SetWavespeed (tel, wavespeed[tent->els[elnr]]);
          int eli = ndomains > 1 ? macroel[tent->els[elnr]] : 0;
          CalcTentElEval (tent->els[elnr], tent, tel, sir, slh,
                          sol.Range (eli * nbasis, (eli + 1) * nbasis),
                          topdshapes[elnr]);
        }
    }); // end loop over tents
    // cout<<"solved from " << timeshift;
    timeshift += tps->GetSlabHeight ();
    // cout<<" to " << timeshift<<endl;
  }

  template <int D>
  template <typename TFUNC>
  void TWaveTents<D>::CalcTentEl (int elnr, const Tent *tent,
                                  ScalarMappedElement<D + 1> &tel,
                                  TFUNC LocalWavespeed,
                                  SIMD_IntegrationRule &sir, LocalHeap &slh,
                                  SliceMatrix<> elmat, FlatVector<> elvec,
                                  SliceMatrix<SIMD<double>> simddshapes)
  {
    static Timer tint1 ("tent top calcshape");
    static Timer tint2 ("tent top AAt");
    static Timer tint3 ("tent top bilinearform");

    HeapReset hr (slh);
    // double wavespeed = tel.GetWavespeed();
    size_t snip = sir.Size () * nsimd;
    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces

    SIMD_STMappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (elnr, slh),
                                                 -1, slh);
    SIMD_MappedIntegrationRule<D, D> smir_fix (sir, ma->GetTrafo (elnr, slh),
                                               slh);
    for (size_t imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ().Range (0, D) = smir_fix[imip].Point ().Range (0, D);

    Mat<D + 1> vert;     // vertices of tent face
    Vec<D + 1> linbasis; // coeffs for linear face fct
    FlatVector<SIMD<double>> mirtimes (sir.Size (), slh);

    /// Integration over bot of tent
    // rows of bdbmat / entries in bdbvec correspond to setting up trial
    // functions
    vert = TentFaceVerts (tent, elnr, -1);
    linbasis = vert.Row (D);
    try
      {
        faceint.Evaluate (sir, linbasis, mirtimes);
      }
    catch (ExceptionNOSIMD const &)
      {
        IntegrationRule ir (eltyp, order * 2);
        FlatVector<double> mirt (sir.Size (),
                                 reinterpret_cast<double *> (&mirtimes (0)));
        faceint.Evaluate (ir, linbasis, mirt);
      }

    for (size_t imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    double area = TentFaceArea (vert);
    Vec<D + 1> n = TentFaceNormal (vert, -1);
    FlatVector<> bdbvec ((D + 1) * snip, slh);
    bdbvec = 0;
    for (size_t imip = 0; imip < snip; imip++)
      {
        double weight = sir[imip / nsimd].Weight ()[imip % nsimd] * area;
        bdbvec (D * snip + imip)
            += n (D) * pow (LocalWavespeed (imip), -2) * weight
               * wavefront (elnr, ((!fosystem) + D) * snip + imip);
        for (int d = 0; d < D; d++)
          {
            bdbvec (d * snip + imip)
                += n (D) * weight
                   * wavefront (elnr, ((!fosystem) + d) * snip + imip);
            bdbvec (d * snip + imip)
                -= n (d) * weight
                   * wavefront (elnr, ((!fosystem) + D) * snip + imip);
            bdbvec (D * snip + imip)
                -= n (d) * weight
                   * wavefront (elnr, ((!fosystem) + d) * snip + imip);
          }
      }
    tel.CalcDShape (smir, simddshapes);
    FlatMatrix<> bbmat (nbasis, (D + 1) * snip,
                        reinterpret_cast<double *> (&simddshapes (0, 0)));
    elvec -= bbmat * bdbvec;

    // stabilization to recover second order solution
    if (!fosystem)
      {
        FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
        tel.CalcShape (smir, simdshapes);
        for (size_t imip = 0; imip < sir.Size (); imip++)
          simdshapes.Col (imip) *= sqrt (area * sir[imip].Weight ());
        AddABt (simdshapes, simdshapes, elmat);
        for (size_t imip = 0; imip < sir.Size (); imip++)
          simdshapes.Col (imip) *= sqrt (area * sir[imip].Weight ());
        FlatMatrix<> shapes (nbasis, snip,
                             reinterpret_cast<double *> (&simdshapes (0, 0)));
        elvec += shapes * wavefront.Row (elnr).Range (0, snip);
      }

    /// Integration over top of tent
    vert = TentFaceVerts (tent, elnr, 1);
    linbasis = vert.Row (D);
    try
      {
        faceint.Evaluate (sir, linbasis, mirtimes);
      }
    catch (ExceptionNOSIMD const &)
      {
        IntegrationRule ir (eltyp, order * 2);
        FlatVector<double> mirt (sir.Size (),
                                 reinterpret_cast<double *> (&mirtimes (0)));
        faceint.Evaluate (ir, linbasis, mirt);
      }
    for (size_t imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    tint1.Start ();
    tel.CalcDShape (smir, simddshapes);
    tint1.Stop ();

    tint2.Start ();
    area = TentFaceArea (vert);
    n = TentFaceNormal (vert, 1);
    FlatMatrix<double> bdbmat ((D + 1) * snip, nbasis, slh);
    bdbmat = 0;
    for (size_t imip = 0; imip < snip; imip++)
      {
        double weight = sir[imip / nsimd].Weight ()[imip % nsimd] * area;
        bdbmat.Row (D * snip + imip) += n (D) * pow (LocalWavespeed (imip), -2)
                                        * weight * bbmat.Col (D * snip + imip);
        for (int d = 0; d < D; d++)
          {
            bdbmat.Row (d * snip + imip)
                += n (D) * weight * bbmat.Col (d * snip + imip);
            bdbmat.Row (d * snip + imip)
                -= n (d) * weight * bbmat.Col (D * snip + imip);
            bdbmat.Row (D * snip + imip)
                -= n (d) * weight * bbmat.Col (d * snip + imip);
          }
      }
    tint2.Stop ();
    tint3.Start ();
    elmat += bbmat * bdbmat;
    tint3.Stop ();

    tint3.AddFlops (bbmat.Height () * bbmat.Width () * bdbmat.Width ());
    // cout << "height " << bbmat.Height() << " width " << bbmat.Width()  <<
    // endl;

    // Vec<D+1> nnn = TentFaceNormal(vert,1);
    // static bool displaytenterror = true;
    // if(L2Norm(nnn.Range(0,D))*(wavespeed) > nnn[D] && displaytenterror){
    // cout << endl << "some tents pitched too high" <<
    // endl;displaytenterror=false;
    //}
    // if(TentFaceArea<D>(vert)<smir_fix[0].GetMeasure()[0]&&
    // tent->nbtime[0]==0)
    //{
    // cout << "tent area problem "<<endl<< TentFaceArea<D>(vert)<< " smaller
    // thamn " <<smir_fix[0].GetMeasure()[0]<< " ratio " <<
    // smir_fix[0].GetJacobiDet()[0]/TentFaceArea<D>(vert) << " vert " << endl;
    // for(int iii=0;iii<D+1;iii++) cout<<vert.Row(iii)<<endl;
    //}
  }

  template <int D>
  void TWaveTents<D>::CalcTentBndEl (int surfel, const Tent *tent,
                                     ScalarMappedElement<D + 1> &tel,
                                     SIMD_IntegrationRule &sir, LocalHeap &slh,
                                     SliceMatrix<> elmat, FlatVector<> elvec)
  {
    HeapReset hr (slh);
    size_t snip = sir.Size () * nsimd;

    // get vertices of tent face
    Mat<D + 1> vert = TentFaceVerts (tent, surfel, 0);
    // build normal vector
    Vec<D + 1> n;
    n = -TentFaceNormal (vert, 0);

    if (D == 1) // D=1 special case
      {
        n[0] = sgn_nozero<int> (tent->vertex - tent->nbv[0]);
        n[D] = 0;
      }
    // build mapping to physical boundary simplex
    Mat<D + 1, D> map;
    for (int i = 0; i < D; i++)
      map.Col (i) = vert.Col (i + 1) - vert.Col (0);
    Vec<D + 1> shift = vert.Col (0);

    SIMD_STMappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (0, slh),
                                                 -1, slh);
    for (size_t imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ()
          = map * sir[imip].operator Vec<D, SIMD<double>> () + shift;

    FlatMatrix<SIMD<double>> simddshapes ((D + 1) * nbasis, sir.Size (), slh);
    tel.CalcDShape (smir, simddshapes);
    FlatMatrix<double> bbmat (
        nbasis, (D + 1) * snip,
        reinterpret_cast<double *> (&simddshapes (0, 0)));

    for (size_t imip = 0; imip < smir.Size (); imip++)
      smir[imip].Point ()[D] += timeshift;
    double area = TentFaceArea (vert);
    FlatMatrix<double> bdbmat ((D + 1) * snip, nbasis, slh);
    bdbmat = 0;
    FlatVector<> bdbvec ((D + 1) * snip, slh);
    bdbvec = 0;
    if (ma->GetMaterial (ElementId (BND, surfel)) == "neumann")
      {
        // had trouble evaluating the normal when passing bndc as n*grad(u),
        // which is why it is now done here. Not useful for actual Neumann
        // bndc, only for testing
        FlatMatrix<SIMD<double>> bdeval (D, sir.Size (), slh);
        bddatum->Evaluate (smir, bdeval);
        double beta = 0.5;
        for (size_t imip = 0; imip < snip; imip++)
          for (int r = 0; r < (D + 1); r++)
            for (int d = 0; d < (D + 1); d++)
              {
                bdbmat.Row (r * snip + imip)
                    += (d < D ? -n (d) * beta : 1.0) * (-n (r))
                       * sir[imip / nsimd].Weight ()[imip % nsimd] * area
                       * bbmat.Col (d * snip + imip);
                bdbvec (d * snip + imip)
                    += (d < D ? -n (d) * beta : -1.0) * (-n (r))
                       * bdeval (r, imip / nsimd)[imip % nsimd]
                       * sir[imip / nsimd].Weight ()[imip % nsimd] * area;
              }
        elmat += bbmat * bdbmat;
        elvec += bbmat * bdbvec;
      }
    else
      { // dirichlet
        FlatMatrix<SIMD<double>> bdeval (1, sir.Size (), slh);
        bddatum->Evaluate (smir, bdeval);
        FlatMatrix<SIMD<double>> wavespeed (1, sir.Size (), slh);
        auto localwavespeedcf = make_shared<ConstantCoefficientFunction> (1)
                                / (this->wavespeedcf * this->wavespeedcf);
        localwavespeedcf->Evaluate (smir, wavespeed);
        double alpha = 0.5;

        for (size_t imip = 0; imip < snip; imip++)
          {
            double weight = area * sir[imip / nsimd].Weight ()[imip % nsimd];
            bdbmat.Row (D * snip + imip)
                += weight * alpha * wavespeed (0, imip / nsimd)[imip % nsimd]
                   * bbmat.Col (D * snip + imip);
            bdbvec (D * snip + imip)
                -= weight * alpha * wavespeed (0, imip / nsimd)[imip % nsimd]
                   * bdeval (0, imip / nsimd)[imip % nsimd];
            for (int d = 0; d < D; d++)
              {
                bdbmat.Row (D * snip + imip)
                    -= n (d) * weight * bbmat.Col (d * snip + imip);
                bdbvec (d * snip + imip)
                    -= n (d) * weight * bdeval (0, imip / nsimd)[imip % nsimd];
              }
          }
        elmat += bbmat * bdbmat;
        elvec -= bbmat * bdbvec;
      }
  }

  template <int D>
  void
  TWaveTents<D>::CalcTentMacroEl (int fnr, const Array<int> &elnums,
                                  std::unordered_map<int, int> &macroel,
                                  const Tent *tent,
                                  ScalarMappedElement<D + 1> &tel,
                                  SIMD_IntegrationRule &sir, LocalHeap &slh,
                                  SliceMatrix<> elmat, FlatVector<>)
  {
    HeapReset hr (slh);
    size_t snip = sir.Size () * nsimd;

    // Array<int> fnums;

    // get vertices of tent face
    // Mat<D+1> vert = TentFaceVerts(tent, surfel, 0);
    // auto sel_verts = ma->GetElVertices(ElementId(BND,elnr));
    Mat<D + 1, D + 1> vert;
    Array<int> sel_verts (D);
    ma->GetFacetPNums (fnr, sel_verts);
    // vert.Col(0).Range(0,D) = ma->GetPoint<D>(tent->vertex);
    vert.Col (0) = ma->GetPoint<D> (tent->vertex);
    vert (D, 0) = tent->tbot;
    for (int n = 0; n < D; n++)
      {
        // vert.Col(n+1).Range(0,D) = ma->GetPoint<D>(sel_verts[n]);
        vert.Col (n + 1) = ma->GetPoint<D> (sel_verts[n]);
        vert (D, n + 1) = tent->vertex == sel_verts[n]
                              ? tent->ttop
                              : tent->nbtime[tent->nbv.Pos (sel_verts[n])];
      }

    // build normal vector
    IntegrationRule ir (eltyp, order * 2);
    MappedIntegrationRule<D, D> mir_fix (ir, ma->GetTrafo (elnums[0], slh),
                                         slh);
    auto fnums = ma->GetElFacets (elnums[0]);
    mir_fix.ComputeNormalsAndMeasure (eltyp, fnums.Pos (fnr));
    auto n = mir_fix[0].GetNV ();

    // build mapping to physical boundary simplex
    Mat<D + 1, D> map;
    for (int i = 0; i < D; i++)
      map.Col (i) = vert.Col (i + 1) - vert.Col (0);
    Vec<D + 1> shift = vert.Col (0);

    SIMD_STMappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (0, slh),
                                                 -1, slh);
    for (size_t imip = 0; imip < snip; imip++)
      smir[imip].Point ()
          = map * sir[imip].operator Vec<D, SIMD<double>> () + shift;

    FlatMatrix<> *bbmat[2];

    SetWavespeed (tel, this->wavespeed[elnums[0]]);
    FlatMatrix<SIMD<double>> simddshapes1 ((D + 1) * nbasis, sir.Size (), slh);
    tel.CalcDShape (smir, simddshapes1);
    bbmat[0]
        = new FlatMatrix<> (nbasis, (D + 1) * snip,
                            reinterpret_cast<double *> (&simddshapes1 (0, 0)));

    SetWavespeed (tel, this->wavespeed[elnums[1]]);
    FlatMatrix<SIMD<double>> simddshapes2 ((D + 1) * nbasis, sir.Size (), slh);
    tel.CalcDShape (smir, simddshapes2);
    bbmat[1]
        = new FlatMatrix<> (nbasis, (D + 1) * snip,
                            reinterpret_cast<double *> (&simddshapes2 (0, 0)));

    FlatMatrix<> *bdbmat[4];
    for (int i = 0; i < 4; i++)
      {
        bdbmat[i] = new FlatMatrix<> ((D + 1) * snip, nbasis, slh);
        *bdbmat[i] = 0;
      }
    // double alpha = 0;
    // double beta = 0;

    double area = TentFaceArea (vert);
    for (size_t imip = 0; imip < snip; imip++)
      {
        double weight = sir[imip / nsimd].Weight ()[imip % nsimd] * area;
        for (int el = 0; el < 4; el++)
          {
            for (int d = 0; d < D; d++)
              {
                bdbmat[el]->Row (d * snip + imip)
                    -= pow (-1, el / 2) * 0.5 * n (d) * weight
                       * bbmat[el % 2]->Col (D * snip + imip);
                bdbmat[el]->Row (D * snip + imip)
                    -= pow (-1, el / 2) * 0.5 * n (d) * weight
                       * bbmat[el % 2]->Col (d * snip + imip);
              }
          }
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
  void TWaveTents<D>::CalcTentElEval (int elnr, const Tent *tent,
                                      ScalarMappedElement<D + 1> &tel,
                                      SIMD_IntegrationRule &sir,
                                      LocalHeap &slh, SliceVector<> sol,
                                      SliceMatrix<SIMD<double>> simddshapes)
  {
    HeapReset hr (slh);
    size_t snip = sir.Size () * nsimd;
    ScalarFE<eltyp, 1> faceint; // linear basis for tent faces

    SIMD_STMappedIntegrationRule<D, D + 1> smir (sir, ma->GetTrafo (elnr, slh),
                                                 -1, slh);
    SIMD_MappedIntegrationRule<D, D> smir_fix (sir, ma->GetTrafo (elnr, slh),
                                               slh);
    for (size_t imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point ().Range (0, D) = smir_fix[imip].Point ().Range (0, D);

    Mat<D + 1> v = TentFaceVerts (tent, elnr, 1);
    Vec<D + 1> linbasis = v.Row (D);
    FlatVector<SIMD<double>> mirtimes (sir.Size (), slh);
    try
      {
        faceint.Evaluate (sir, linbasis, mirtimes);
      }
    catch (ExceptionNOSIMD const &)
      {
        IntegrationRule ir (eltyp, order * 2);
        FlatVector<double> mirt (sir.Size (),
                                 reinterpret_cast<double *> (&mirtimes (0)));
        faceint.Evaluate (ir, linbasis, mirt);
      }
    for (size_t imip = 0; imip < sir.Size (); imip++)
      smir[imip].Point () (D) = mirtimes[imip];

    FlatMatrix<SIMD<double>> simdshapes (nbasis, sir.Size (), slh);
    // FlatMatrix<SIMD<double>> simddshapes((D+1)*nbasis,sir.Size(),slh);
    if (!fosystem)
      tel.CalcShape (smir, simdshapes);
    // tel.CalcDShape(smir,simddshapes);
    FlatMatrix<> dshapes (nbasis, (D + 1) * snip,
                          reinterpret_cast<double *> (&simddshapes (0, 0)));
    FlatMatrix<> shapes (nbasis, snip,
                         reinterpret_cast<double *> (&simdshapes (0, 0)));
    if (!fosystem)
      wavefront.Row (elnr).Range (0, snip) = Trans (shapes) * sol;
    wavefront.Row (elnr).Range (snip * (!fosystem),
                                snip * (!fosystem) + snip * (D + 1))
        = Trans (dshapes) * sol;
  }

  // returns matrix where cols correspond to vertex coordinates of the
  // space-time element
  template <int D>
  Mat<D + 1, D + 1>
  TWaveTents<D>::TentFaceVerts (const Tent *tent, int elnr, int top)
  {
    Mat<D + 1, D + 1> v;
    if (top == 0) // boundary element
      {
        auto sel_verts = ma->GetElVertices (ElementId (BND, elnr));
        // v.Col(0).Range(0,D) = ma->GetPoint<D>(tent->vertex);
        v.Col (0) = ma->GetPoint<D> (tent->vertex);
        v (D, 0) = tent->tbot;
        for (int n = 0; n < D; n++)
          {
            // v.Col(n+1).Range(0,D) = ma->GetPoint<D>(sel_verts[n]);
            v.Col (n + 1) = ma->GetPoint<D> (sel_verts[n]);
            v (D, n + 1) = tent->vertex == sel_verts[n]
                               ? tent->ttop
                               : tent->nbtime[tent->nbv.Pos (sel_verts[n])];
          }
      }
    else // top or bot of tent
      {
        IVec<D + 1> vnr = ma->GetElVertices (ElementId (VOL, elnr));
        // determine linear basis function coeffs to use for tent face
        for (size_t ivert = 0; ivert < vnr.Size (); ivert++)
          {
            // v.Col(ivert).Range(0,D) = ma->GetPoint<D>(vnr[ivert]);
            v.Col (ivert) = ma->GetPoint<D> (vnr[ivert]);
            if (vnr[ivert] == tent->vertex)
              v (D, ivert) = top == 1 ? tent->ttop : tent->tbot;
            else
              for (size_t k = 0; k < tent->nbv.Size (); k++)
                if (vnr[ivert] == tent->nbv[k])
                  v (D, ivert) = tent->nbtime[k];
          }
      }
    return v;
  }

  template <int D> double TWaveTents<D>::TentFaceArea (Mat<D + 1, D + 1> ve)
  {
    // note: we scale by factor 2 or 6 (D=2,3) the measure of the reference
    // element. Already computed in GetWeight/GetMeasure of any MappedIP
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
  Vec<D + 1> TWaveTents<D>::TentFaceNormal (Mat<D + 1, D + 1> v, int top)
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
          Vec<3> ba = v.Col (0) - v.Col (1);
          Vec<3> baa = v.Col (0) - v.Col (2);
          normv = Cross (ba, baa);
          break;
        }
      case 3:
        {
          for (int d = 1; d < D + 1; d++)
            v.Col (d) = v.Col (0) - v.Col (d);

          for (int i = 0; i < D + 1; i++)
            {
              Mat<D, D> pS;
              for (int k = 0, c = 0; k < D + 1; k++)
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
  Matrix<> TWaveTents<D>::MakeWavefront (shared_ptr<CoefficientFunction> cf,
                                         double time)
  {
    LocalHeap lh (1000 * 1000 * 1000, "make wavefront", 1);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    size_t snip = sir.Size () * nsimd;
    Matrix<> wf (ma->GetNE (), snip * cf->Dimension ());
    for (size_t elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_STMappedIntegrationRule<D, D + 1> smir (
            sir, ma->GetTrafo (elnr, lh), -1, lh);
        SIMD_MappedIntegrationRule<D, D> smir_fix (
            sir, ma->GetTrafo (elnr, lh), lh);
        for (size_t imip = 0; imip < sir.Size (); imip++)
          {
            smir[imip].Point ().Range (0, D)
                = smir_fix[imip].Point ().Range (0, D);
            smir[imip].Point ()[D] = time;
          }
        FlatMatrix<SIMD<double>> bdeval (cf->Dimension (), smir.Size (), lh);
        bdeval = 0;
        cf->Evaluate (smir, bdeval);
        for (size_t imip = 0; imip < snip; imip++)
          for (size_t d = 0; d < cf->Dimension (); d++)
            wf (elnr, d * snip + imip)
                = bdeval (d, imip / nsimd)[imip % nsimd];
      }
    return wf;
  }

  template <int D>
  double TWaveTents<D>::Error (Matrix<> wavefront, Matrix<> wavefront_corr)
  {
    LocalHeap lh (1000 * 1000 * 1000, "error", 1);
    double error = 0;
    SIMD_IntegrationRule sir (eltyp, order * 2);
    size_t snip = sir.Size () * nsimd;
    for (size_t elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);

        for (size_t imip = 0; imip < snip; imip++)
          {
            IntegrationRule ir (eltyp, 0);
            MappedIntegrationPoint<D, D> mip (
                ir[0], ma->GetTrafo (ElementId (0), lh));
            for (int d = 0; d < D; d++)
              mip.Point ()[d] = smir[imip / nsimd].Point ()[d][0];
            // mip.Point()[D] = timeshift;
            for (int d = 0; d < D + 1; d++)
              error += pow (wavespeedcf->Evaluate (mip), -2 * (d == D))
                       * pow (wavefront (elnr, ((!fosystem) + d) * snip + imip)
                                  - wavefront_corr (
                                      elnr, ((!fosystem) + d) * snip + imip),
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
  double TWaveTents<D>::L2Error (Matrix<> wavefront, Matrix<> wavefront_corr)
  {
    LocalHeap lh (1000 * 1000 * 1000, "l2error", 1);
    double l2error = 0;
    SIMD_IntegrationRule sir (eltyp, order * 2);
    size_t snip = sir.Size () * nsimd;
    for (size_t elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);
        for (size_t imip = 0; imip < snip; imip++)
          {
            l2error += (wavefront (elnr, imip) - wavefront_corr (elnr, imip))
                       * (wavefront (elnr, imip) - wavefront_corr (elnr, imip))
                       * smir[imip / nsimd].GetWeight ()[imip % nsimd];
          }
      }
    return sqrt (l2error);
  }

  template <int D> double TWaveTents<D>::Energy (Matrix<> wavefront)
  {
    double energy = 0;
    LocalHeap lh (1000 * 1000 * 1000, "energy", 1);
    SIMD_IntegrationRule sir (eltyp, order * 2);
    size_t snip = sir.Size () * nsimd;
    for (size_t elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        HeapReset hr (lh);
        SIMD_MappedIntegrationRule<D, D> smir (sir, ma->GetTrafo (elnr, lh),
                                               lh);
        FlatMatrix<SIMD<double>> wavespeed (1, smir.Size (), lh);
        wavespeedcf->Evaluate (smir, wavespeed);
        for (size_t imip = 0; imip < snip; imip++)
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
  void TWaveTents<D>::SwapIfGreater (T &a, T &b)
  {
    if (a < b)
      {
        T tmp (a);
        a = b;
        b = tmp;
      }
  }

  template <int D> double TWaveTents<D>::TentAdiam (const Tent *tent)
  {
    LocalHeap lh (1000 * 1000 * 1000);
    int vnumber = tent->nbv.Size ();
    // double c = wavespeed[tent->els[0]];
    // for(auto el : tent->els) c = max(c,wavespeed[el]);
    IntegrationRule ir (ma->GetElType (ElementId (VOL, 0)), 0);
    ElementTransformation &trafo = ma->GetTrafo (0, lh);
    MappedIntegrationPoint<D, D> mip (ir[0], trafo);
    // Vec<D> v1 = ma->GetPoint<D>(tent->vertex);
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

  template <int D> double TWaveTents<D>::MaxAdiam ()
  {
    double h = 0.0;
    RunParallelDependency (tps->tent_dependency, [&] (int tentnr) {
      const Tent *tent = &tps->GetTent (tentnr);
      h = max (h, TentAdiam (tent));
    });
    return h;
  }

  template <int D>
  inline int TWaveTents<D>::MakeMacroEl (const Array<int> &tentel,
                                         std::unordered_map<int, int> &macroel)
  {
    // TODO fix if macro elements do not share faces
    int nrmacroel = 0;
    for (size_t i = 0; i < tentel.Size (); i++)
      {
        size_t j = 0;
        while (wavespeed[tentel[i]] != wavespeed[tentel[j]])
          j++;
        if (j == i)
          macroel[tentel[i]] = nrmacroel++;
        else
          macroel[tentel[i]] = macroel[tentel[j]];
      }
    return nrmacroel;
  }

  template <int D>
  void TWaveTents<D>::GetFacetSurfaceElement (shared_ptr<MeshAccess> ma,
                                              int fnr, Array<int> &selnums)
  {
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
          for (size_t i : Range (ma->GetNSE ()))
            {
              auto sel = ElementId (BND, i); //<-- ist das Oberflächenelement
              auto fnums = ma->GetElFacets (sel);
              if (fnr == fnums[0]) //<-- das ist dann die Facet-Nr., also unter
                                   // allen Facets im Mesh, kannst du dir wo
                                   // speichern
                {
                  selnums.Append (i);
                  break;
                }
            }
        }
      }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <int D> void QTWaveTents<D>::Propagate ()
  {
    LocalHeap lh (1000 * 1000 * 1000, "QT tents", 1);

    shared_ptr<MeshAccess> ma = this->ma;
    SIMD_IntegrationRule sir (eltyp, this->order * 2);

    // cout << "solving qt " << (this->tps)->GetNTents() << " tents in " << D
    // << "+1 dimensions..." << endl;

    RunParallelDependency ((this->tps)->tent_dependency, [&] (int tentnr) {
      LocalHeap slh = lh.Split (); // split to threads
      const Tent *tent = &(this->tps)->GetTent (tentnr);

      Vec<D + 1> center;
      center.Range (0, D) = ma->GetPoint<D> (tent->vertex);
      center[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
      double tentsize = TentXdiam (tent);

      // QTWaveFE<D> tel(GGder, BBder, this->order, center, tentsize);
      CSR basismat = this->basis.Basis (this->order, center, tentsize);
      int nbasis = this->nbasis;
      ScalarMappedElement<D + 1> tel (nbasis, this->order, basismat, ET_TET,
                                      center, 1.0 / tentsize);

      FlatMatrix<> elmat (nbasis, slh);
      FlatVector<> elvec (nbasis, slh);
      elmat = 0;
      elvec = 0;

      Array<FlatMatrix<SIMD<double>>> topdshapes (tent->els.Size ());
      for (auto &tds : topdshapes)
        tds.AssignMemory ((D + 1) * nbasis, sir.Size (), slh);

      // Integrate top and bottom space-like tent faces
      for (size_t elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          HeapReset hr (slh);
          SIMD_MappedIntegrationRule<D, D> smir_fix (
              sir, ma->GetTrafo (tent->els[elnr], slh), slh);
          FlatMatrix<SIMD<double>> lwavespeed (1, sir.Size (), slh);
          (this->wavespeedcf)->Evaluate (smir_fix, lwavespeed);

          this->CalcTentEl (
              tent->els[elnr], tent, tel,
              [&] (int imip) {
                return lwavespeed (0, imip / nsimd)[imip % nsimd];
              },
              sir, slh, elmat, elvec, topdshapes[elnr]);
        }

      for (auto fnr : tent->internal_facets)
        {
          Array<int> elnums;
          ma->GetFacetElements (fnr, elnums);
          Array<int> selnums;
          if (elnums.Size () == 1)
            this->GetFacetSurfaceElement (ma, fnr, selnums);

          // Integrate boundary tent
          if (elnums.Size () == 1 && selnums.Size () == 1)
            this->CalcTentBndEl (selnums[0], tent, tel, sir, slh, elmat,
                                 elvec);
        }

      // integrate volume of tent here
      for (size_t elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          /// Integration over bot and top volume of tent element
          for (int part = -1; part <= 1; part += 2)
            {
              // if(tent->nbtime[elnr]==(part<0?tent->tbot:tent->ttop))
              // continue;
              HeapReset hr (slh);
              const ELEMENT_TYPE vol_eltyp = (D == 2) ? ET_TET : ET_TRIG;
              SIMD_IntegrationRule vsir (vol_eltyp, this->order * 2);

              Vec<D + 1> shift;
              shift.Range (0, D) = ma->GetPoint<D> (tent->vertex);
              shift[D] = (tent->ttop - tent->tbot) / 2 + tent->tbot;
              // shift[D] = tent->nbtime[elnr];
              // shift[D] = part==1 ? tent->nbtime[elnr] : tent->tbot;
              Mat<D + 1> map
                  = this->TentFaceVerts (tent, tent->els[elnr], part);

              for (int i = 0; i < D + 1; i++)
                map.Col (i) -= shift;
              double vol = abs (
                  Det (map)); // no need for this * (D==2?6.0:2.0); bc of Det
              if (vol < 10e-16)
                continue;

              SIMD_MappedIntegrationRule<D + 1, D + 1> vsmir (
                  vsir, ma->GetTrafo (tent->els[elnr], slh), -1, slh);
              for (size_t imip = 0; imip < vsir.Size (); imip++)
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
              tel.CalcDDWaveOperator (vsmir, simdddshapes, wavespeed);
              for (size_t imip = 0; imip < vsir.Size (); imip++)
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
              for (size_t j = 0; j < vsir.Size ();
                   j++) // auto ws : wavespeed.AsVector())
                for (size_t i = 0; i < nsimd; i++)
                  cmax = max (cmax, abs (wavespeed (0, j)[i]));
              FlatMatrix<SIMD<double>> mu (1, vsir.Size (), slh);
              localwavespeedcf
                  = make_shared<ConstantCoefficientFunction> (tentsize / cmax)
                    * this->wavespeedcf * this->wavespeedcf;
              localwavespeedcf->Evaluate (vsmir, mu);
              FlatMatrix<SIMD<double>> simdddshapescor ((D + 1) * nbasis,
                                                        vsir.Size (), slh);
              tel.CalcDDWaveOperator (vsmir, simdddshapescor, wavespeed, mu);
              AddABt (FlatMatrix<SIMD<double>> (nbasis, (D + 1) * vsir.Size (),
                                                &simdddshapescor (0, 0)),
                      simdddshapes2, elmat);
            }
        }

      // solve
      Solve (elmat, elvec);
      FlatVector<> sol (nbasis, &elvec (0));

      // eval solution on top of tent
      for (size_t elnr = 0; elnr < tent->els.Size (); elnr++)
        {
          this->CalcTentElEval (tent->els[elnr], tent, tel, sir, slh, sol,
                                topdshapes[elnr]);
        }
    }); // end loop over tents
    // cout<<"solved from " << this->timeshift;
    this->timeshift += (this->tps)->GetSlabHeight ();
    // cout<<" to " << this->timeshift<<endl;
  }

  template <int D> double QTWaveTents<D>::TentXdiam (const Tent *tent)
  {
    int vnumber = tent->nbv.Size ();
    double diam = 0;

    shared_ptr<MeshAccess> ma = this->ma;
    for (int k = 0; k < vnumber; k++)
      {
        Vec<D> v1 = ma->GetPoint<D> (tent->nbv[k]);
        Vec<D> v2 = ma->GetPoint<D> (tent->vertex);
        diam = max (diam, sqrt (L2Norm2 (v1 - v2)));
        for (int j = k; j < vnumber; j++)
          {
            v2 = ma->GetPoint<D> (tent->nbv[j]);
            diam = max (diam, sqrt (L2Norm2 (v1 - v2)));
          }
      }

    return diam;
  }

  template class TWaveTents<1>;
  template class TWaveTents<2>;
  template class TWaveTents<3>;

  template class QTWaveTents<1>;
  template class QTWaveTents<2>;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>

template <class T, int D>
void DeclareETClass (py::module &m, std::string typestr)
{
  using PyETclass = T;
  std::string pyclass_name = typestr;
  py::class_<PyETclass, shared_ptr<PyETclass>, TrefftzTents> (
      m, pyclass_name.c_str ())
      //.def(py::init<>())
      .def ("MakeWavefront", &PyETclass::MakeWavefront)
      .def ("GetWavefront", &PyETclass::GetWavefront)
      .def ("Error", &PyETclass::Error)
      .def ("L2Error", &PyETclass::L2Error)
      .def ("Energy", &PyETclass::Energy)
      .def ("MaxAdiam", &PyETclass::MaxAdiam)
      .def ("LocalDofs", &PyETclass::LocalDofs)
      .def ("GetOrder", &PyETclass::GetOrder)
      .def ("GetSpaceDim", &PyETclass::GetSpaceDim)
      .def ("GetInitmesh", &PyETclass::GetInitmesh);
}

void ExportTWaveTents (py::module m)
{
  py::class_<TrefftzTents, shared_ptr<TrefftzTents>> (m, "TrefftzTents")
      .def ("Propagate", &TrefftzTents::Propagate, "Solve tent slab")
      .def ("SetInitial", &TrefftzTents::SetInitial, "Set initial condition")
      .def ("SetBoundaryCF", &TrefftzTents::SetBoundaryCF,
            "Set boundary condition");

  DeclareETClass<TWaveTents<1>, 1> (m, "TWaveTents1");
  DeclareETClass<TWaveTents<2>, 2> (m, "TWaveTents2");
  DeclareETClass<TWaveTents<3>, 3> (m, "TWaveTents3");
  DeclareETClass<QTWaveTents<1>, 1> (m, "QTWaveTents1");
  DeclareETClass<QTWaveTents<2>, 2> (m, "QTWaveTents2");

  m.def (
      "TWave",
      [] (int order, shared_ptr<TentPitchedSlab> tps,
          shared_ptr<CoefficientFunction> wavespeedcf,
          shared_ptr<CoefficientFunction> BBcf) -> shared_ptr<TrefftzTents> {
        shared_ptr<TrefftzTents> tr;
        int D = (tps->ma)->GetDimension ();
        if (!BBcf)
          {
            if (D == 1)
              tr = make_shared<TWaveTents<1>> (order, tps, wavespeedcf);
            else if (D == 2)
              tr = make_shared<TWaveTents<2>> (order, tps, wavespeedcf);
            else if (D == 3)
              tr = make_shared<TWaveTents<3>> (order, tps, wavespeedcf);
          }
        else
          {
            if (D == 1)
              tr = make_shared<QTWaveTents<1>> (order, tps, wavespeedcf, BBcf);
            else if (D == 2)
              tr = make_shared<QTWaveTents<2>> (order, tps, wavespeedcf, BBcf);
          }
        return tr;
      },
      R"mydelimiter(
                Create solver for acoustiv wave equation on tent-pitched mesh.

                :param order: Polynomial order of the Trefftz space.
                :param tps: Tent-pitched slab.
                :param wavespeedcf: PDE Coefficient
                :param BB: PDE Coefficient
            )mydelimiter",
      py::arg ("order"), py::arg ("tps"), py::arg ("wavespeedcf"),
      py::arg ("BBcf") = nullptr);
}
#endif // NGS_PYTHON
