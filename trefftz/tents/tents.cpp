#include "tents.hpp"
#ifdef USE_NETGEN_GUI
#include "myvisual.hpp"
#endif

constexpr ELEMENT_TYPE ET_ (int D)
{
  return (D == 1) ? ET_SEGM : (D == 2) ? ET_TRIG : ET_TET;
}

// dirty hack since we done have the numproc
int additional_intorder = 1;

template <int DIM>
void TentPitchedSlab<DIM>::PitchTents (double dt, double wavespeed,
                                       LocalHeap &lh)
{
  auto cf = make_shared<ConstantCoefficientFunction> (wavespeed);
  PitchTents (dt, cf, lh);
}
template <int DIM>
void TentPitchedSlab<DIM>::PitchTents (
    double dt, shared_ptr<CoefficientFunction> wavespeed, LocalHeap &lh)
{
  cout << "pitch tents ... ";
  Timer t ("pitch tents");
  RegionTimer reg (t);
  // Timer tminlevel("pitch tents - find minlevel");
  // Timer tloop1("pitch tents - loop1");
  // Timer tloop2("pitch tents - loop2");
  // Timer tloop3("pitch tents - loop3");
  // Timer tcalc("calc matrices");

  // element-wise maximal wave-speeds
  Array<double> cmax (ma->GetNE ());
  // cmax = wavespeed;

  // compute edge-based max time-differences
  Array<double> edge_refdt (ma->GetNEdges ());
  edge_refdt = 1e99; // max-double ??
  // check if edge is contained in mesh
  BitArray fine_edges (ma->GetNEdges ());
  fine_edges.Clear ();
  for (Ngs_Element el : ma->Elements (VOL))
    {
      ElementId ei = ElementId (el);
      ELEMENT_TYPE eltype = ma->GetElType (ei);
      IntegrationRule ir (eltype, 0);
      ElementTransformation &trafo = ma->GetTrafo (ei, lh);
      MappedIntegrationPoint<DIM, DIM> mip (ir[0], trafo);
      cmax[el.Nr ()] = wavespeed->Evaluate (mip) + 2.0;

      for (int e : el.Edges ())
        {
          auto pnts = ma->GetEdgePNums (e);
          auto v1 = pnts[0], v2 = pnts[1];
          double len = L2Norm (ma->template GetPoint<DIM> (v1)
                               - ma->template GetPoint<DIM> (v2));
          edge_refdt[e] = min (edge_refdt[e], len / cmax[el.Nr ()]);
          fine_edges.Set (e);
        }
    }

  Array<double> vertex_refdt (ma->GetNV ());
  vertex_refdt = 1e99;
  for (int e : IntRange (0, ma->GetNEdges ()))
    if (fine_edges.Test (e))
      {
        auto vts = ma->GetEdgePNums (e);
        int v1 = vts[0], v2 = vts[1];
        vertex_refdt[v1] = min (vertex_refdt[v1], edge_refdt[e]);
        vertex_refdt[v2] = min (vertex_refdt[v2], edge_refdt[e]);
      }

  Array<double> tau (ma->GetNV ()); // advancing front
  tau = 0.0;
  Array<double> ktilde (ma->GetNV ());
  ktilde = vertex_refdt;

  // 'list' of ready vertices
  Array<int> ready_vertices;
  Array<bool> vertex_ready (ma->GetNV ()); // quick check if vertex is ready
  for (int i = 0; i < ma->GetNV (); i++)
    ready_vertices.Append (i);
  vertex_ready = true;

  Array<bool> complete_vertices (ma->GetNV ());
  complete_vertices = false;

  // build vertex2edge and vertex2vertex tables

  TableCreator<int> create_v2e, create_v2v;
  for (; !create_v2e.Done (); create_v2e++, create_v2v++)
    {
      for (int e : IntRange (0, ma->GetNEdges ()))
        if (fine_edges.Test (e))
          {
            auto vts = ma->GetEdgePNums (e);
            int v1 = vts[0], v2 = vts[1];
            create_v2v.Add (v1, v2);
            create_v2v.Add (v2, v1);

            create_v2e.Add (v1, e);
            create_v2e.Add (v2, e);
          }
    }

  Table<int> v2v = create_v2v.MoveTable ();
  Table<int> v2e = create_v2e.MoveTable ();

  Array<int> latest_tent (ma->GetNV ()), vertices_level (ma->GetNV ());
  latest_tent = -1;
  vertices_level = 0;

  while (ready_vertices.Size ())
    {
      int minlevel = 1000;
      int posmin = 0;
      // tminlevel.Start();
      for (int i = 0; i < ready_vertices.Size (); i++)
        if (vertices_level[ready_vertices[i]] < minlevel)
          {
            minlevel = vertices_level[ready_vertices[i]];
            posmin = i;
          }
      // tminlevel.Stop();
      // tloop1.Start();
      int vi = ready_vertices[posmin];
      ready_vertices.DeleteElement (posmin);
      vertex_ready[vi] = false;

      // cout << "\rconstruct tent at vertex " << vi << setw(5) << flush;

      // advance by ktilde:

      Tent *tent = new Tent;
      tent->vertex = vi;
      tent->tbot = tau[vi];
      tent->ttop = min (dt, tau[vi] + ktilde[vi]);
      tent->level = vertices_level[vi]; // 0;
      tau[vi] = tent->ttop;

      for (int nb : v2v[vi])
        {
          tent->nbv.Append (nb);
          tent->nbtime.Append (tau[nb]);
          if (vertices_level[nb] < tent->level + 1)
            vertices_level[nb] = tent->level + 1;
          if (latest_tent[nb] != -1)
            {
              tents[latest_tent[nb]]->dependent_tents.Append (
                  tents.Size ()); // my tent number
              // tent->level = max2 (tent->level,
              // tents[latest_tent[nb]]->level+1);
            }
        }
      latest_tent[vi] = tents.Size ();
      vertices_level[vi]++;
      // tloop1.Stop();
      // tloop2.Start();
      if (DIM == 1)
        tent->edges.Append (
            vi); // vertex itself represents the only internal edge/facet
      else if (DIM == 2)
        for (int e : v2e[vi])
          tent->edges.Append (e);
      else
        {
          ArrayMem<int, 4> fpnts;
          for (auto elnr : ma->GetVertexElements (vi))
            for (auto f : ma->GetElement (ElementId (VOL, elnr)).Faces ())
              {
                ma->GetFacetPNums (f, fpnts);
                if (fpnts.Contains (vi) && !tent->edges.Contains (f))
                  tent->edges.Append (f);
              }
        }

      /*
      for (int f = 0; f < ma->GetNFacets(); f++)
        {
          ArrayMem<int,4> fpnts;
          ma->GetFacetPNums(f, fpnts);
          if (fpnts.Contains(vi))
            tent->edges.Append(f);
        }
      */

      ma->GetVertexElements (vi, tent->els);
      // tloop2.Stop();
      // tloop3.Start();
      // update max step for neighbours:
      for (int nb : v2v[vi])
        {
          if (tau[nb] >= dt)
            continue;

          double kt = 1e99;
          for (int nb2_index : v2v[nb].Range ())
            {
              int nb2 = v2v[nb][nb2_index];
              double kt1 = tau[nb2] - tau[nb] + edge_refdt[v2e[nb][nb2_index]];
              kt = min (kt, kt1);
              // if (tau[nb]+kt > dt) kt = dt-tau[nb];
            }

          ktilde[nb] = kt;
          if (kt > 0.5 * vertex_refdt[nb]) //  || (tau[nb]+kt > dt) )
                                           // if (!ready_vertices.Contains(nb))
            if (!vertex_ready[nb])
              {
                ready_vertices.Append (nb);
                vertex_ready[nb] = true;
              }
        }
      tents.Append (tent);
      // if(DIM == 2){
      // 	netgen::AddUserVisualizationObject (new MyDraw(*ma, tents));
      // 	Ng_Redraw();
      // 	getchar();
      // }
      // tloop3.Stop();
    }
  // cout << endl;
  // build dependency graph
  TableCreator<int> create_dag (tents.Size ());
  for (; !create_dag.Done (); create_dag++)
    {
      for (int i : tents.Range ())
        for (int d : tents[i]->dependent_tents)
          create_dag.Add (i, d);
    }
  tent_dependency = create_dag.MoveTable ();
  cout << "done " << endl;
}

template <int DIM>
double
TentPitchedSlab<DIM>::GetTentHeight (int vertex, Array<int> &els,
                                     FlatArray<int> nbv, Array<double> &tau,
                                     Array<double> &cmax, LocalHeap &lh)
{
  HeapReset hr (lh);

  ma->GetVertexElements (vertex, els);
  // cout << "tau[vertex] = " <<  tau[vertex] << endl;
  double height = 1e99;
  for (int j : Range (els.Size ()))
    {
      ElementId ej (VOL, els[j]);
      auto vnums = ma->GetElVertices (ej);
      ElementTransformation &trafo = ma->GetTrafo (ej, lh);

      //       const DGFiniteElement<DIM> & fel =
      // 	static_cast<const DGFiniteElement<DIM>&> (fes->FESpace::GetFE
      // (ej, lh));

      ELEMENT_TYPE eltype = ma->GetElType (ej);

      ScalarFE<ET_ (DIM), 1> fe_nodal;
      Vector<> shape_nodal (fe_nodal.GetNDof ());
      Matrix<> dshape_nodal (fe_nodal.GetNDof (), DIM);

      Vector<> coef_top (fe_nodal.GetNDof ()); // coefficient of tau_top(X)
      // coef_top.SetSize(fe_nodal.GetNDof());

      int pos = -1;
      for (int k = 0; k < vnums.Size (); k++)
        {
          if (vnums[k] == vertex) // central vertex
            {
              coef_top (k) = 0.0;
              pos = k;
            }
          else
            {
              for (int l = 0; l < nbv.Size (); l++)
                if (nbv[l] == vnums[k])
                  coef_top (k) = tau[nbv[l]];
            }
        }
      // cout << "pos = " << pos << endl;
      // cout << "coef_top = " << endl << coef_top << endl;
      IntegrationRule ir (eltype, 1); // 2*fel.Order());

      MappedIntegrationPoint<DIM, DIM> mip (ir[0],
                                            trafo); // for non-curved elements
      fe_nodal.CalcMappedDShape (mip, dshape_nodal);
      // solve quadratic (in)equality
      // cout << "grad = " << dshape_nodal << endl;
      Matrix<> gradphisqr = dshape_nodal * Trans (dshape_nodal);
      // cout << "mat = " << endl << gradphisqr << endl;
      Vector<> eigenvals (fe_nodal.GetNDof ());
      LapackEigenValuesSymmetric (gradphisqr, eigenvals);
      // cout << "eigenvals = " << eigenvals << endl;
      double alpha
          = InnerProduct (dshape_nodal.Row (pos), dshape_nodal.Row (pos));
      Vec<DIM> temp = Trans (dshape_nodal) * coef_top;
      double beta = 2.0 * InnerProduct (temp, dshape_nodal.Row (pos));
      double gamma = 0.0;
      for (int k = 0; k < fe_nodal.GetNDof (); k++)
        gamma += coef_top (k) * InnerProduct (temp, dshape_nodal.Row (k));
      gamma -= 1.0 / (cmax[ej.Nr ()] * cmax[ej.Nr ()]);

      double temp2 = sqrt (beta * beta - 4.0 * alpha * gamma);
      double sol1 = (-beta + temp2) / (2.0 * alpha);
      double sol2 = (-beta - temp2) / (2.0 * alpha);
      // cout << "alpha = " << alpha << endl;
      // cout << "beta = " << beta << endl;
      // cout << "gamma = " << gamma << endl;
      // cout << "tmax = " << sol1 << ", " << sol2 << endl;
      if (sol1 > sol2)
        {
          // cout << "old = " << tau[vertex] << endl;
          // cout << "new = " << sol1 << endl;
          height = min (height, sol1);
        }
      else
        {
          // cout << "old = " << tau[vertex] << endl;
          // cout << "new = " << sol2 << endl;
          height = min (height, sol2);
        }
      // for(int l = 0; l < coef_top.Size(); l++)
      // if(l!=pos) cout << "tau[nb] = " << coef_top(l) << endl;
      // coef_top(pos) = sol1;
      // temp = Trans(dshape_nodal) * coef_top;
      // cout << "|gradphi|^2 = " << InnerProduct(temp,temp) << ", ";
      // coef_top(pos) = sol2;
      // temp = Trans(dshape_nodal) * coef_top;
      // cout << InnerProduct(temp,temp) << endl;
    }
  height -= tau[vertex];
  // cout << "height = " << height << endl;
  // if(height < 0.0) cout << "height = " << height << " negative
  // !!!!!!!!!!!!!!!!!!!!!!!" << endl;
  return height;
}

template <int DIM>
void TentPitchedSlab<DIM>::PitchTentsGradient (double dt, double wavespeed,
                                               LocalHeap &lh)
{
  cout << "pitch tents ... ";
  Timer t ("pitch tents");
  RegionTimer reg (t);
  // Timer tminlevel("pitch tents - find minlevel");
  // Timer tloop1("pitch tents - loop1");
  // Timer tloop2("pitch tents - loop2");
  // Timer tloop3("pitch tents - loop3");
  // Timer tcalc("calc matrices");

  // element-wise maximal wave-speeds
  Array<double> cmax (ma->GetNE ());
  cmax = wavespeed;

  // compute edge-based max time-differences
  Array<double> edge_refdt (ma->GetNEdges ());
  edge_refdt = 1e99; // max-double ??
  // check if edge is contained in mesh
  BitArray fine_edges (ma->GetNEdges ());
  fine_edges.Clear ();
  for (Ngs_Element el : ma->Elements (VOL))
    {
      for (int e : el.Edges ())
        {
          int v1, v2;
          // ma->GetEdgePNums (e, v1, v2);
          std::tie (v1, v2) = ma->GetEdgePNums (e);
          double len = L2Norm (ma->template GetPoint<DIM> (v1)
                               - ma->template GetPoint<DIM> (v2));
          edge_refdt[e] = min (edge_refdt[e], len / cmax[el.Nr ()]);
          fine_edges.Set (e);
        }
    }

  Array<double> vertex_refdt (ma->GetNV ());
  vertex_refdt = 1e99;
  for (int e : IntRange (0, ma->GetNEdges ()))
    if (fine_edges.Test (e))
      {
        auto pnts = ma->GetEdgePNums (e);
        auto v1 = pnts[0], v2 = pnts[1];
        vertex_refdt[v1] = min (vertex_refdt[v1], edge_refdt[e]);
        vertex_refdt[v2] = min (vertex_refdt[v2], edge_refdt[e]);
      }

  *testout << "old: " << endl << vertex_refdt << endl;
  // build vertex2edge and vertex2vertex tables

  TableCreator<int> create_v2e, create_v2v;
  for (; !create_v2e.Done (); create_v2e++, create_v2v++)
    {
      for (int e : IntRange (0, ma->GetNEdges ()))
        if (fine_edges.Test (e))
          {
            int v1, v2;
            std::tie (v1, v2) = ma->GetEdgePNums (e);
            create_v2v.Add (v1, v2);
            create_v2v.Add (v2, v1);

            create_v2e.Add (v1, e);
            create_v2e.Add (v2, e);
          }
    }

  Table<int> v2v = create_v2v.MoveTable ();
  Table<int> v2e = create_v2e.MoveTable ();

  ArrayMem<int, 30> vels;
  vels.SetSize (0);

  Array<double> tau (ma->GetNV ()); // advancing front
  tau = 0.0;
  Array<double> ktilde (ma->GetNV ());

  // 'list' of ready vertices
  Array<int> ready_vertices;
  Array<bool> vertex_ready (ma->GetNV ()); // quick check if vertex is ready
  for (int i = 0; i < ma->GetNV (); i++)
    {
      ready_vertices.Append (i);
      vertex_refdt[i] = GetTentHeight (i, vels, v2v[i], tau, cmax, lh);
    }
  *testout << "new: " << endl << vertex_refdt << endl;
  ktilde = vertex_refdt;
  vertex_ready = true;

  Array<bool> complete_vertices (ma->GetNV ());
  complete_vertices = false;

  Array<int> latest_tent (ma->GetNV ()), vertices_level (ma->GetNV ());
  latest_tent = -1;
  vertices_level = 0;

  while (ready_vertices.Size ())
    {
      // cout << ready_vertices << endl;
      // cout << tau << endl;
      // cout << complete_vertices << endl;
      // getchar();
      int minlevel = 1000;
      int posmin = 0;
      // tminlevel.Start();
      for (int i = 0; i < ready_vertices.Size (); i++)
        if (vertices_level[ready_vertices[i]] < minlevel)
          {
            minlevel = vertices_level[ready_vertices[i]];
            posmin = i;
          }
      // tminlevel.Stop();
      // tloop1.Start();
      int vi = ready_vertices[posmin];
      ready_vertices.DeleteElement (posmin);
      vertex_ready[vi] = false;

      // cout << "\rconstruct tent at vertex " << vi << setw(5) << flush;

      // advance by ktilde:
      double newktilde = GetTentHeight (vi, vels, v2v[vi], tau, cmax, lh);
      if (newktilde < ktilde[vi])
        {
          cout << "updated" << endl;
          ktilde[vi] = newktilde;
        }
      if (ktilde[vi] < 0.0)
        {
          cout << "vertex " << vi << " : ";
          cout << ktilde[vi] << " ktilde negative !!!!!!!!!!!!!" << endl;
          ktilde[vi] = GetTentHeight (vi, vels, v2v[vi], tau, cmax, lh);
          cout << "new ktilde = " << ktilde[vi] << endl;
          // getchar();
        }

      if (tau[vi] >= dt)
        continue;
      // cout << "dt - tau[vi] = " << dt - tau[vi] << ", ktilde = " <<
      // ktilde[vi] << endl;
      if (ktilde[vi] < 0.5 * vertex_refdt[vi])
        {
          cout << "ktilde[vi] = " << ktilde[vi]
               << ", ref = " << vertex_refdt[vi] << endl;
          continue;
        }
      Tent *tent = new Tent;
      tent->vertex = vi;
      tent->tbot = tau[vi];
      if (tau[vi] + ktilde[vi] > dt)
        {
          tent->ttop = dt;
          complete_vertices[vi] = true;
        }
      else
        {
          // cout << "max diff = " << dt*1e-8 << endl;
          if (dt - tau[vi] - ktilde[vi] < dt * 1e-8)
            tent->ttop = tau[vi] + 0.9 * ktilde[vi];
          else
            tent->ttop = tau[vi] + ktilde[vi];
        }
      tent->level = vertices_level[vi]; // 0;
      tau[vi] = tent->ttop;

      for (int nb : v2v[vi])
        {
          tent->nbv.Append (nb);
          tent->nbtime.Append (tau[nb]);
          if (vertices_level[nb] < tent->level + 1)
            vertices_level[nb] = tent->level + 1;
          if (latest_tent[nb] != -1)
            {
              tents[latest_tent[nb]]->dependent_tents.Append (
                  tents.Size ()); // my tent number
              // tent->level = max2 (tent->level,
              // tents[latest_tent[nb]]->level+1);
            }
        }
      latest_tent[vi] = tents.Size ();
      vertices_level[vi]++;
      // tloop1.Stop();
      // tloop2.Start();
      if (DIM == 1)
        tent->edges.Append (
            vi); // vertex itself represents the only internal edge/facet
      else if (DIM == 2)
        for (int e : v2e[vi])
          tent->edges.Append (e);
      else
        {
          ArrayMem<int, 4> fpnts;
          for (auto elnr : ma->GetVertexElements (vi))
            for (auto f : ma->GetElement (ElementId (VOL, elnr)).Faces ())
              {
                ma->GetFacetPNums (f, fpnts);
                if (fpnts.Contains (vi) && !tent->edges.Contains (f))
                  tent->edges.Append (f);
              }
        }

      /*
      for (int f = 0; f < ma->GetNFacets(); f++)
        {
          ArrayMem<int,4> fpnts;
          ma->GetFacetPNums(f, fpnts);
          if (fpnts.Contains(vi))
            tent->edges.Append(f);
        }
      */

      ma->GetVertexElements (vi, tent->els);
      // tloop2.Stop();
      // tloop3.Start();
      // cout << "old ktilde = " << ktilde[vi] << endl;
      // double maxheightnew = GetTentHeight(*tent,cmax,fes,lh);
      // cout << "new ktilde = " << maxheightnew << endl;
      // update max step for neighbours:
      for (int nb : v2v[vi])
        {
          if (tau[nb] >= dt)
            continue;

          // double kt = 1e99;
          // for (int nb2_index : v2v[nb].Range())
          //   {
          //     int nb2 = v2v[nb][nb2_index];
          //     double kt1 = tau[nb2]-tau[nb]+edge_refdt[v2e[nb][nb2_index]];
          //     kt = min (kt, kt1);
          //     // if (tau[nb]+kt > dt) kt = dt-tau[nb];
          //   }
          //
          // ktilde[nb] = kt;
          // bool nb_ready = false;
          // for (int nb2 : v2v[nb])
          //   if(vertex_ready[nb2]) nb_ready = true;
          // if(nb_ready) continue;

          double temp = ktilde[nb];
          ktilde[nb] = GetTentHeight (nb, vels, v2v[nb], tau, cmax, lh);
          // if(ktilde[nb] - temp < 0)
          //   {
          //     cout << "new step smaller !!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
          //     cout << ktilde[nb] << " < " << temp << endl;
          //     // getchar();
          //   }
          if (ktilde[nb] > 0.5 * vertex_refdt[nb]) //  || (tau[nb]+kt > dt) )
            // if (!ready_vertices.Contains(nb))
            if (!vertex_ready[nb] && !complete_vertices[nb])
              {
                ready_vertices.Append (nb);
                vertex_ready[nb] = true;
                // *testout << "vertex " << nb << " ready to pitch tau = " <<
                // ktilde[nb] << endl;
              }
        }
      tents.Append (tent);
      // if(DIM == 2){
      // 	netgen::AddUserVisualizationObject (new MyDraw(*ma, tents));
      // 	Ng_Redraw();
      // 	getchar();
      // }
      // tloop3.Stop();
    }
  // cout << endl;
  // build dependency graph
  TableCreator<int> create_dag (tents.Size ());
  for (; !create_dag.Done (); create_dag++)
    {
      for (int i : tents.Range ())
        for (int d : tents[i]->dependent_tents)
          create_dag.Add (i, d);
    }
  tent_dependency = create_dag.MoveTable ();
  cout << "done " << endl;
}

template <int DIM>
void TentPitchedSlab<DIM>::SetupTents (
    const shared_ptr<L2HighOrderFESpace> fes, LocalHeap &lh)
{
  spacetime_dofs = 0;
  int cnt_tents = 0;
  cout << "calc data for tents ... ";
  // for (int i = 0; i < tents.Size(); i++)
  ParallelFor (Range (tents), [&] (int i) {
    Array<int> dnums;
    Tent &tent = *tents[i];

    for (int j = 0; j < tent.els.Size (); j++)
      {
        fes->GetDofNrs (tent.els[j], dnums);
        tent.ranges.Append (IntRange (dnums.Size ()) + tent.dofs.Size ());
        tent.dofs += dnums;
        tent.nd_T.Append (dnums.Size ());
      }
    tent.nd = tent.dofs.Size ();
    spacetime_dofs += tent.dofs.Size ();

    TableCreator<int> elfnums_creator (tent.els.Size ());

    for (; !elfnums_creator.Done (); elfnums_creator++)
      {
        for (int j : Range (tent.els))
          {
            auto fnums = ma->GetElFacets (tent.els[j]);
            for (int fnum : fnums)
              if (tent.edges.Pos (fnum) != -1)
                elfnums_creator.Add (j, fnum);
          }
      }
    tent.elfnums = elfnums_creator.MoveTable ();
  });
  /*
  ParallelFor
    (Range(tents), [&] (int i)
     {
// #pragma omp atomic
//        cnt_tents++;
//
// #pragma omp critical
//        cout << "\rcalc data for tent " << cnt_tents << "/" << tents.Size()
<< flush;

       LocalHeap & clh = lh, lh = clh.Split();

       Tent & tent = *tents[i];
       int nels = tent.els.Size();
       tent.gradphi_bot.SetSize(nels);
       tent.gradphi_top.SetSize(nels);
       // tent.agradphi_bot.SetSize(nels);
       // tent.agradphi_top.SetSize(nels);
       tent.delta.SetSize(nels);
       // tent.adelta.SetSize(nels);
       tent.graddelta.SetSize(nels);

       FlatArray<int> elfacets(nels,lh);
       for (int j : Range(nels))
         elfacets[j] = tent.elfnums[j].Size();

       tent.gradphi_facet_bot = Table< Matrix<> >(elfacets);
       tent.gradphi_facet_top = Table< Matrix<> >(elfacets);
       // tent.agradphi_facet_bot = Table< AFlatMatrix<> >(elfacets);
       // tent.agradphi_facet_top = Table< AFlatMatrix<> >(elfacets);
       tent.delta_facet = Table< Vector<double> >(elfacets);
       // tent.adelta_facet = Table< AVector<double> >(elfacets);

       for (int j : Range(nels)) // loop over elts in a tent
         {
           ElementId ej (VOL, tent.els[j]);

           Array<int> dnums;

           fes->GetDofNrs (ej, dnums);
           auto vnums = ma->GetElVertices (ej);

           ElementTransformation & trafo = ma->GetTrafo (ej, lh);

           const DGFiniteElement<DIM> & fel =
             static_cast<const DGFiniteElement<DIM>&> (fes->GetFE (ej, lh));

           ELEMENT_TYPE eltype = fel.ElementType();

           // Set top and bottom tent surfaces using pw.linear FE:

           // ScalarFE< (DIM==2) ? ET_TRIG : ET_SEGM,1> fe_nodal; // for 2d or
1d ScalarFE<ET_(DIM),1> fe_nodal; Vector<> shape_nodal(fe_nodal.GetNDof());
           Matrix<> dshape_nodal(fe_nodal.GetNDof(),DIM);

           Vector<> coef_bot, coef_top; // coefficient of tau_bot(X) and
tau_top(X) coef_bot.SetSize(fe_nodal.GetNDof());
           coef_top.SetSize(fe_nodal.GetNDof());
           for (int k = 0; k < vnums.Size(); k++)
             {
               if (vnums[k] == tent.vertex)  // central vertex
                 {
                   coef_bot(k) = tent.tbot;
                   coef_top(k) = tent.ttop;
                 }
               else
                 for (int l = 0; l < tent.nbv.Size(); l++)
                   if (tent.nbv[l] == vnums[k])
                     coef_bot(k) = coef_top(k) = tent.nbtime[l];
             }

           IntegrationRule ir (eltype, 2*fel.Order());

           tent.gradphi_bot[j].SetSize(ir.Size(),DIM);
           tent.gradphi_top[j].SetSize(ir.Size(),DIM);

           // void * ptr;
           // posix_memalign (&ptr, 64, 8*DIM*(ir.Size()+8));
           // new (&tent.agradphi_top[j]) AFlatMatrix<>(DIM, ir.Size(),
(double*)ptr);
           // posix_memalign (&ptr, 64, 8*DIM*(ir.Size()+8));
           // new (&tent.agradphi_bot[j]) AFlatMatrix<>(DIM, ir.Size(),
(double*)ptr);

           tent.delta[j].SetSize(ir.Size());
           // tent.adelta[j] = AVector<>(ir.Size());
           tent.graddelta[j].SetSize(DIM);

           for (int k = 0; k < ir.Size(); k++)
             {
               MappedIntegrationPoint<DIM,DIM> mip(ir[k], trafo);

               fe_nodal.CalcShape (ir[k], shape_nodal);
               fe_nodal.CalcMappedDShape (mip, dshape_nodal);

               tent.delta[j][k] = InnerProduct(coef_top-coef_bot, shape_nodal);
               tent.gradphi_bot[j].Row(k) = Trans (dshape_nodal) * coef_bot;
               tent.gradphi_top[j].Row(k) = Trans (dshape_nodal) * coef_top;
               if(k==0)
                 tent.graddelta[j] = Trans (dshape_nodal) *
(coef_top-coef_bot);
             }
           // tent.adelta[j] = 0.0;
           // tent.adelta[j] = tent.delta[j];

           // for (int k = 0; k < tent.agradphi_bot[j].VWidth(); k++)
           //   for (int l = 0; l < DIM; l++)
           //     {
           //       tent.agradphi_bot[j].Get(l,k) = 0.0;
           //       tent.agradphi_top[j].Get(l,k) = 0.0;
           //     }

           // tent.agradphi_bot[j] = 0.0;
           // tent.agradphi_top[j] = 0.0;
           // for (int k = 0; k < tent.agradphi_bot[j].Width(); k++)
           //   for (int l = 0; l < DIM; l++)
           //     {
           //       tent.agradphi_bot[j](l,k) = tent.gradphi_bot[j](k,l);
           //       tent.agradphi_top[j](l,k) = tent.gradphi_top[j](k,l);
           //     }

           for (int k : tent.elfnums[j].Range())
             {
               const DGFiniteElement<DIM-1> & felfacet =
                 static_cast<const DGFiniteElement<DIM-1>&> (fes->GetFacetFE
(tent.elfnums[j][k], lh)); IntegrationRule ir(felfacet.ElementType(),
2*felfacet.Order()+additional_intorder);

               tent.gradphi_facet_bot[j][k].SetSize(ir.Size(), DIM);
               tent.gradphi_facet_top[j][k].SetSize(ir.Size(), DIM);
               tent.delta_facet[j][k].SetSize(ir.Size());
               // tent.adelta_facet[j][k] = AVector<>(ir.Size());

               auto vnums = ma->GetElVertices (ej);
               Facet2ElementTrafo transform(fel.ElementType(), vnums);

               int loc_facetnr = 0; // local facet number of element j
               // if( facetdata[tent.elfnums[j][k]]->elnr[0] == ej.Nr() )
               // 	 loc_facetnr =
facetdata[tent.elfnums[j][k]]->facetnr[0];
               // else
               //	 loc_facetnr =
facetdata[tent.elfnums[j][k]]->facetnr[1];

               auto fnums = ma->GetElFacets(ej);
               for (int l = 0; l < fnums.Size(); l++)
                 if (fnums[l] == tent.elfnums[j][k]) loc_facetnr = l;

               for (int l = 0; l < ir.Size(); l++)
                 {
                   MappedIntegrationPoint<DIM,DIM> mip(transform(loc_facetnr,
ir[l]), trafo);

                   fe_nodal.CalcShape (mip.IP(), shape_nodal);
                   fe_nodal.CalcMappedDShape (mip, dshape_nodal);

                   tent.delta_facet[j][k][l] = InnerProduct(coef_top-coef_bot,
shape_nodal); tent.gradphi_facet_bot[j][k].Row(l) = Trans (dshape_nodal) *
coef_bot; tent.gradphi_facet_top[j][k].Row(l) = Trans (dshape_nodal) *
coef_top;
                 }
               // tent.adelta_facet[j][k] = 0.0;
               // tent.adelta_facet[j][k] = tent.delta_facet[j][k];
               //
               // posix_memalign (&ptr, 64, 8*DIM*(ir.Size()+8));
               // new (&tent.agradphi_facet_top[j][k]) AFlatMatrix<>(DIM,
ir.Size(), (double*)ptr);
               // posix_memalign (&ptr, 64, 8*DIM*(ir.Size()+8));
               // new (&tent.agradphi_facet_bot[j][k]) AFlatMatrix<>(DIM,
ir.Size(), (double*)ptr);
               //
               // tent.agradphi_facet_bot[j][k] = 0.0;
               // tent.agradphi_facet_top[j][k] = 0.0;
               // for (int l = 0; l < tent.agradphi_facet_bot[j][k].Width();
l++)
               //   for (int m = 0; m < DIM; m++)
               //     {
               //       tent.agradphi_facet_bot[j][k](m,l) =
tent.gradphi_facet_bot[j][k](l,m);
               //       tent.agradphi_facet_top[j][k](m,l) =
tent.gradphi_facet_top[j][k](l,m);
               //     }
             }
         } // loop over elements of tent
     }); // loop over all tents
  */
  cout << "done" << endl;
  cout << "number of tents = " << tents.Size () << endl;
}

TempTentData::TempTentData (const Tent &tent, const FESpace &fes,
                            const MeshAccess &ma, LocalHeap &lh)
    : fei (tent.els.Size (), lh), iri (tent.els.Size (), lh),
      miri (tent.els.Size (), lh), trafoi (tent.els.Size (), lh),
      mesh_size (tent.els.Size (), lh), agradphi_bot (tent.els.Size (), lh),
      agradphi_top (tent.els.Size (), lh), adelta (tent.els.Size (), lh),
      felpos (tent.edges.Size (), lh), firi (tent.edges.Size (), lh),
      mfiri1 (tent.edges.Size (), lh), mfiri2 (tent.edges.Size (), lh),
      agradphi_botf1 (tent.edges.Size (), lh),
      agradphi_topf1 (tent.edges.Size (), lh),
      agradphi_botf2 (tent.edges.Size (), lh),
      agradphi_topf2 (tent.edges.Size (), lh),
      anormals (tent.edges.Size (), lh), adelta_facet (tent.edges.Size (), lh)
{
  int dim = ma.GetDimension ();
  FlatArray<BaseScalarFiniteElement *> fe_nodal (tent.els.Size (), lh);
  FlatArray<FlatVector<double>> coef_delta (tent.els.Size (), lh);
  FlatArray<FlatVector<double>> coef_top (tent.els.Size (), lh);
  FlatArray<FlatVector<double>> coef_bot (tent.els.Size (), lh);

  for (int i = 0; i < tent.els.Size (); i++)
    {
      ElementId ei (VOL, tent.els[i]);
      fei[i] = &fes.GetFE (ei, lh);
      iri[i] = new (lh)
          SIMD_IntegrationRule (fei[i]->ElementType (), 2 * fei[i]->Order ());
      trafoi[i] = &ma.GetTrafo (ei, lh);
      miri[i] = &(*trafoi[i]) (*iri[i], lh);

      mesh_size[i] = pow (fabs ((*miri[i])[0].GetJacobiDet ()[0]),
                          1.0 / miri[i]->DimElement ());

      int nipt = miri[i]->Size ();
      agradphi_bot[i].AssignMemory (dim, nipt, lh);
      agradphi_top[i].AssignMemory (dim, nipt, lh);
      adelta[i].AssignMemory (nipt, lh);

      // BaseScalarFiniteElement * fe_nodal;
      switch (dim)
        {
        case 1:
          fe_nodal[i] = new (lh) ScalarFE<ET_SEGM, 1> ();
          break;
        case 2:
          fe_nodal[i] = new (lh) ScalarFE<ET_TRIG, 1> ();
          break;
        default:
          fe_nodal[i] = new (lh) ScalarFE<ET_TET, 1> ();
        }

      // VectorMem<4,double> coef_bot(fe_nodal[i]->GetNDof()),
      // coef_top(fe_nodal[i]->GetNDof());
      coef_top[i].AssignMemory (fe_nodal[i]->GetNDof (), lh);
      coef_bot[i].AssignMemory (fe_nodal[i]->GetNDof (), lh);
      auto vnums = ma.GetElVertices (ei);
      for (int k = 0; k < vnums.Size (); k++)
        {
          auto pos = tent.nbv.Pos (vnums[k]);
          if (pos != tent.nbv.ILLEGAL_POSITION)
            coef_bot[i](k) = coef_top[i](k) = tent.nbtime[pos];
          else
            {
              coef_bot[i](k) = tent.tbot;
              coef_top[i](k) = tent.ttop;
            }
        }
      fe_nodal[i]->EvaluateGrad (*miri[i], coef_top[i], agradphi_top[i]);
      fe_nodal[i]->EvaluateGrad (*miri[i], coef_bot[i], agradphi_bot[i]);

      coef_delta[i].AssignMemory (fe_nodal[i]->GetNDof (), lh);
      coef_delta[i] = coef_top[i] - coef_bot[i];
      fe_nodal[i]->Evaluate (*iri[i], coef_delta[i], adelta[i]);
    }

  for (int i = 0; i < tent.edges.Size (); i++)
    {
      int order = 0;
      INT<2> loc_facetnr;

      ArrayMem<int, 2> elnums;
      ma.GetFacetElements (tent.edges[i], elnums);

      SIMD_IntegrationRule *simd_ir_facet;

      felpos[i] = INT<2, size_t> (size_t (-1));
      for (int j : Range (elnums.Size ()))
        {
          felpos[i][j] = tent.els.Pos (elnums[j]);
          if (felpos[i][j] != size_t (-1))
            {
              int elorder = fei[felpos[i][j]]->Order ();
              if (elorder > order)
                order = elorder;

              auto fnums = ma.GetElFacets (elnums[j]);
              for (int k : Range (fnums.Size ()))
                if (fnums[k] == tent.edges[i])
                  loc_facetnr[j] = k;

              auto &trafo = *trafoi[felpos[i][j]];

              auto vnums = ma.GetElVertices (elnums[j]);
              Facet2ElementTrafo transform (trafo.GetElementType (), vnums);

              auto etfacet = ElementTopology::GetFacetType (
                  trafo.GetElementType (), loc_facetnr[j]);
              if (j == 0)
                {
                  simd_ir_facet = new (lh) SIMD_IntegrationRule (
                      etfacet, 2 * order + additional_intorder);
                  simd_ir_facet->SetIRX (
                      nullptr); // quick fix to avoid usage of TP elements
                                // (slows down)
                }

              // auto & simd_ir_facet_vol = transform(loc_facetnr[j],
              // *simd_ir_facet, lh);
              firi[i][j] = &transform (loc_facetnr[j], *simd_ir_facet, lh);
              if (j == 0)
                {
                  mfiri1[i] = &trafo (*firi[i][j], lh);
                  mfiri1[i]->ComputeNormalsAndMeasure (trafo.GetElementType (),
                                                       loc_facetnr[j]);

                  anormals[i].AssignMemory (dim, mfiri1[i]->Size (), lh);
                  adelta_facet[i].AssignMemory (mfiri1[i]->Size (), lh);
                  agradphi_botf1[i].AssignMemory (dim, mfiri1[i]->Size (), lh);
                  agradphi_topf1[i].AssignMemory (dim, mfiri1[i]->Size (), lh);

                  for (size_t k : Range (mfiri1[i]->Size ()))
                    {
                      switch (dim)
                        {
                        case 1:
                          anormals[i].Col (k) = static_cast<const SIMD<
                              DimMappedIntegrationPoint<1>> &> (
                                                    (*mfiri1[i])[k])
                                                    .GetNV ();
                          break;
                        case 2:
                          anormals[i].Col (k) = static_cast<const SIMD<
                              DimMappedIntegrationPoint<2>> &> (
                                                    (*mfiri1[i])[k])
                                                    .GetNV ();
                          break;
                        default:
                          anormals[i].Col (k) = static_cast<const SIMD<
                              DimMappedIntegrationPoint<3>> &> (
                                                    (*mfiri1[i])[k])
                                                    .GetNV ();
                        }
                    }
                  size_t elpos = felpos[i][j];
                  fe_nodal[elpos]->Evaluate (*firi[i][j], coef_delta[elpos],
                                             adelta_facet[i]);
                  fe_nodal[elpos]->EvaluateGrad (*mfiri1[i], coef_bot[elpos],
                                                 agradphi_botf1[i]);
                  fe_nodal[elpos]->EvaluateGrad (*mfiri1[i], coef_top[elpos],
                                                 agradphi_topf1[i]);
                }
              else
                {
                  mfiri2[i] = &trafo (*firi[i][j], lh);
                  mfiri2[i]->ComputeNormalsAndMeasure (trafo.GetElementType (),
                                                       loc_facetnr[j]);
                  agradphi_botf2[i].AssignMemory (dim, mfiri2[i]->Size (), lh);
                  agradphi_topf2[i].AssignMemory (dim, mfiri2[i]->Size (), lh);
                  size_t elpos = felpos[i][j];
                  fe_nodal[elpos]->EvaluateGrad (*mfiri2[i], coef_bot[elpos],
                                                 agradphi_botf2[i]);
                  fe_nodal[elpos]->EvaluateGrad (*mfiri2[i], coef_top[elpos],
                                                 agradphi_topf2[i]);
                }
            }
        }
    }
}

template <int DIM>
void TentPitchedSlab<DIM>::CheckTents (
    const shared_ptr<L2HighOrderFESpace> fes, LocalHeap &lh)
{
  cout << "call CheckTents" << endl;
  double maxgrad = 0.0;
  ParallelFor (Range (tents), [&] (int i) {
    LocalHeap &clh = lh, lh = clh.Split ();
    Tent &tent = *tents[i];
    *testout << "tent " << i << ":" << endl;
    for (int j : Range (tent.els.Size ())) // loop over elts in a tent
      {
        // *testout << "el " << j << ":" << endl;
        int height = tent.gradphi_top[j].Height ();
        auto grad = tent.gradphi_top[j].Row (0);
        auto norm = L2Norm (grad);
        *testout << norm << endl;
        if (norm > maxgrad)
          maxgrad = norm;
      }
  });
  cout << "maximal gradient = " << maxgrad << endl;
}

template <int DIM> void TentPitchedSlab<DIM>::VTKOutputTents (string filename)
{
  ;
}

template <int DIM>
void TentPitchedSlab<DIM>::GetTentData (Array<int> &tentdata,
                                        Array<double> &tentbot,
                                        Array<double> &tenttop, int &nlevels)
{
  nlevels = 0;
  tentdata.SetAllocSize (4 * tents.Size ());
  tenttop.SetAllocSize (4 * tents.Size ());
  tentbot.SetAllocSize (4 * tents.Size ());

  for (int i : Range (tents))
    {
      Tent &tent = *tents[i];
      for (int el : Range (tent.els))
        {
          tentdata.Append (i);
          tentdata.Append (tent.level);
          tentdata.Append (tent.vertex);
          tentdata.Append (tent.els[el]);
          if (tent.level > nlevels)
            nlevels = tent.level;

          auto verts = ma->GetElVertices (ElementId (VOL, tent.els[el]));
          for (auto v : verts)
            {
              if (v == tent.vertex)
                tentbot.Append (tent.tbot);
              else
                {
                  auto pos = tent.nbv.Pos (v);
                  tentbot.Append (tent.nbtime[pos]);
                }
            }
          tentbot.Append (tent.ttop);
          tenttop.Append ({ tent.ttop, 0.0, 0.0, 0.0 });
        }
    }
  nlevels += 1;
}

void VTKOutputTents (shared_ptr<MeshAccess> maptr, Array<Tent *> &tents,
                     string filename)
{
  const MeshAccess &ma = *maptr;
  ofstream out (filename + ".vtk");
  Array<Vec<3>> points;
  Array<INT<4>> cells;
  Array<int> level, tentnr;
  int ptcnt = 0;

  for (int i : Range (tents))
    {
      int firstpt = ptcnt;
      Tent &tent = *tents[i];
      Vec<2> pxy = ma.GetPoint<2> (tent.vertex);
      points.Append (Vec<3> (pxy (0), pxy (1), tent.tbot));
      points.Append (Vec<3> (pxy (0), pxy (1), tent.ttop));
      INT<4> tet (ptcnt, ptcnt + 1, 0, 0);
      ptcnt += 2;

      for (int elnr : tent.els)
        {
          Ngs_Element el = ma.GetElement (ElementId (VOL, elnr));

          for (int v : el.Vertices ())
            if (v != tent.vertex)
              {
                pxy = ma.GetPoint<2> (v);
                points.Append (
                    Vec<3> (pxy (0), pxy (1), tent.nbtime[tent.nbv.Pos (v)]));
              }
          for (int j = 2; j < 4; j++)
            tet[j] = ptcnt++;

          cells.Append (tet);
        }
      for (int j = firstpt; j < ptcnt; j++)
        {
          level.Append (tent.level);
          tentnr.Append (i);
        }
    }

  // header
  out << "# vtk DataFile Version 3.0" << endl;
  out << "vtk output" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;

  out << "POINTS " << points.Size () << " float" << endl;
  for (auto p : points)
    out << p << endl;

  out << "CELLS " << cells.Size () << " " << 5 * cells.Size () << endl;
  for (auto c : cells)
    out << 4 << " " << c << endl;

  out << "CELL_TYPES " << cells.Size () << endl;
  for (auto c : cells)
    out << "10 " << endl;

  out << "CELL_DATA " << cells.Size () << endl;
  out << "POINT_DATA " << points.Size () << endl;

  out << "FIELD FieldData " << 2 << endl;

  out << "tentlevel"
      << " 1 " << level.Size () << " float" << endl;
  for (auto i : level)
    out << i << " ";
  out << endl;

  out << "tentnumber"
      << " 1 " << tentnr.Size () << " float" << endl;
  for (auto i : tentnr)
    out << i << " ";
  out << endl;
}

template class TentPitchedSlab<1>;
template class TentPitchedSlab<2>;
template class TentPitchedSlab<3>;
