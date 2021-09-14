#include "tents.hpp"
#include <limits>
#include <h1lofe.hpp> // seems needed for ScalarFE (post 2021-06-22 NGSolve update)



///////////////////// GradPhiCoefficientFunction ///////////////////////////

void GradPhiCoefficientFunction::GenerateCode(Code &code, FlatArray<int> inputs, int index) const
{
  auto dims = Dimensions();

  string header = "\n\
    {flatmatrix} {values};\n\
    ProxyUserData * {ud} = (ProxyUserData*)mir.GetTransformation().userdata;\n\
    {\n\
      // class GradPhiCoefficientFunction\n ";
  header+=
    "      if ({ud}->fel) {\n";
  if(code.is_simd) {
    header += "auto x = {ud}->GetAMemory ({this});\n";
    header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
  } else {
    header += "auto x = {ud}->GetMemory ({this});\n";
    header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
  }
  header +=  "}\n";
  header += "}\n";

  string body = "";
  TraverseDimensions( dims, [&](int ind, int i, int j) {
      body += Var(index, i,j).Declare("{scal_type}", 0.0);
      string values = "{values}";
      if(code.is_simd)
	values += "(" + ToLiteral(ind) + ",i)";
      else
	values += "(i," + ToLiteral(ind) + ")";
      body += Var(index, i,j).Assign(CodeExpr(values), false);
    });

  std::map<string,string> variables;
  variables["ud"] = "tmp_"+ToLiteral(index)+"_0";
  variables["this"] = "reinterpret_cast<CoefficientFunction*>("+code.AddPointer(this)+")";

  variables["flatmatrix"] = code.is_simd ? "FlatMatrix<SIMD<double>>" : "FlatMatrix<double>";
  variables["values"] = Var("values", index).S();

  string scal_type = "double";
  if(code.is_simd) scal_type = "SIMD<" + scal_type + ">";
  variables["scal_type"] = scal_type;

  code.header += Code::Map(header, variables);
  code.body += Code::Map(body, variables);
}

/////////////////// Tent meshing ///////////////////////////////////////////

constexpr ELEMENT_TYPE EL_TYPE(int DIM)
{
  return DIM == 1 ? ET_SEGM : DIM == 2 ? ET_TRIG : ET_TET;
}//this assumes that there is only one type of element per mesh

template <int DIM>
bool TentPitchedSlab::PitchTents(const double dt, const bool calc_local_ct, const double global_ct)
{
  if(has_been_pitched)
    {
      tents.DeleteAll();
    }
  if(cmax == nullptr)
    {
      throw std::logic_error("Wavespeed has not been set!");
    }
  this->dt = dt; // set it so that GetSlabHeight can return it
  TentSlabPitcher * slabpitcher = [this]() ->TentSlabPitcher* {
    switch (this->method)
      {
      case ngstents::EVolGrad:
        return new VolumeGradientPitcher<DIM>(this->ma, vmap);
        break;
      case ngstents::EEdgeGrad:
        return new EdgeGradientPitcher<DIM>(this->ma, vmap);
      default:
        cout << "Trying to pitch tent without setting a pitching method." << endl;
        return nullptr;
        break;
      }
  }();
  if(!slabpitcher) return false;
  cout << "Created slab pitcher"<<endl;
  //calc wavespeed for each element and perhaps other stuff (i..e, calculating edge gradients, checking fine edges, etc)
  Table<int> v2v, v2e;
  std::tie(v2v,v2e) = slabpitcher->InitializeMeshData<DIM>(lh,cmax, calc_local_ct, global_ct);
  cout << "Initialised mesh data" << endl;
  
  Array<double> tau(ma->GetNV());  // advancing front values at vertices
  tau = 0.0;

  
  slabpitcher->ComputeVerticesReferenceHeight(v2v, v2e, tau, lh);
  cout << "Computed reference heights" << endl;
  // max time increase allowed at vertex, depends on tau of neighbors
  //at the beginning the advancing front is at a constant t=0
  //so ktilde can be set as vertex_refdt
  Array<double> ktilde = slabpitcher->GetVerticesReferenceHeight();

  //array containing the latest tent in which the vertex was included
  Array<int> latest_tent(ma->GetNV());
  //array containing the current level of the tent to be pitched at vertex vi
  Array<int> vertices_level(ma->GetNV());
  latest_tent = -1;
  //every vertex start at level 0
  vertices_level = 0;

  //for an advance to be considered good, dt >= factor * refdt
  double adv_factor{0.5};
  //whether to reset the adv_factor to its initial value after populating ready_vertices
  constexpr bool reset_adv_factor = true;
  // array of vertices ready for pitching a tent
  Array<int> ready_vertices;
  // array for checking if a given vertex is ready
  BitArray vertex_ready(ma->GetNV());
  vertex_ready.Clear();
  bool slab_complete{false};
  //array for checking if a given vertex is complete (tau[vi] = dt)
  BitArray complete_vertices(ma->GetNV());
  complete_vertices.Clear();
  //numerical tolerance
  const double num_tol = std::numeric_limits<double>::epsilon() * dt;

  while ( !slab_complete )
    {
      cout << "Setting ready vertices" << endl;
      const bool found_vertices =
        slabpitcher->GetReadyVertices(adv_factor,reset_adv_factor,ktilde,complete_vertices, vertex_ready,ready_vertices);
      //no possible vertex in which a tent could be pitched was found
      if(!found_vertices) break;
      // ---------------------------------------------
      // Main loop: constructs one tent each iteration
      // ---------------------------------------------
      cout << "Pitching tents..." << endl;      
      while (ready_vertices.Size())
        {
          int minlevel, posmin;
          std::tie(minlevel,posmin) =
            slabpitcher->PickNextVertexForPitching(ready_vertices,ktilde,vertices_level);
          nlayers = max(minlevel,nlayers);
          //vertex index at which the current tent is being pitched
          const int vi = ready_vertices[posmin];
          ready_vertices.DeleteElement(posmin);
          vertex_ready.Clear(vi);

          //current tent
          Tent * tent = new Tent(vmap);
          tent->vertex = vi;
          tent->tbot = tau[vi];

          const auto new_ttop = tau[vi] + ktilde[vi];
          if(dt - new_ttop > num_tol)
            {//not close to the end of the time slab
              tent->ttop = new_ttop;
            }
          else
            {//vertex is complete
              tent->ttop = dt;
              complete_vertices.SetBit(vi);
            }
          //let us ignore this for now
          // else if(new_ttop >= dt)
          //   {//vertex is complete
          //     tent->ttop = dt;
          //     complete_vertices[vi] = true;
          //   }
          // else
          //   {//vertex is really close to the end of time slab.
          //     //in this scenario, we might want to pitch a lower
          //     //tent to avoid numerical issues with degenerate tents
          //     tent->ttop = ktilde[vi] * 0.75 + tau[vi];
          //   }
          
          tent->level = vertices_level[vi]; // 0;
          tau[vi] = tent->ttop;
          ktilde[vi] = 0;//assuming that ktilde[vi] was the maximum advance

          //add neighboring vertices and update their level
          for (int nb : v2v[vi])
            {
              nb = vmap[nb]; // only update main vertex if periodic
              tent->nbv.Append (nb);
              tent->nbtime.Append (tau[nb]);
              //update level of vertices if needed
              if(vertices_level[nb] < tent->level + 1)
                vertices_level[nb] = tent->level + 1;
              // tent number is just array index in tents
              if (latest_tent[nb] != -1)
                tents[latest_tent[nb]]->dependent_tents.Append (tents.Size());
            }
          latest_tent[vi] = tents.Size();
          vertices_level[vi]++;

          // Set tent internal facets
          if(DIM==1)
            // vertex itself represents the only internal edge/facet
            tent->internal_facets.Append (vi);
          else if (DIM == 2)
            for (int e : v2e[vi]) tent->internal_facets.Append (e);
          else
            {
              // DIM == 3 => internal facets are faces
              //points contained in a given facet
              ArrayMem<int,4> fpnts;
              ArrayMem<int,30> vertex_els;
              slabpitcher->GetVertexElements(vi,vertex_els);
              for (auto elnr : vertex_els)
                for (auto f : ma->GetElement(ElementId(VOL,elnr)).Faces())
                  {
                    //get facet vertices
                    ma->GetFacetPNums(f, fpnts);
                    for (auto f_v : fpnts)
                      {
                        if (vmap[f_v]  == vi &&
                            !tent->internal_facets.Contains(f))
                          {
                            tent->internal_facets.Append(f);
                            break;
                          }
                      }
                  }
            }          
          slabpitcher->GetVertexElements(vi,tent->els);
          slabpitcher->UpdateNeighbours(vi,adv_factor,v2v,v2e,tau,complete_vertices,
                                        ktilde,vertex_ready,ready_vertices,lh);
          tents.Append (tent);
        }
      //check if slab is complete
      slab_complete = true;
      for(int i = 0; i < ma->GetNV(); i++)
        if(vmap[i] == i)
          if(complete_vertices[i] == false)
            {
              slab_complete = false;
              break;
            }
    }

 
  if(!slab_complete)
    {
      const auto &vrefdt = slabpitcher->GetVerticesReferenceHeight();
      cout << "Error: the algorithm could not pitch the whole slab" << endl;
      int iv;
      for(iv = 0; iv < ma->GetNV(); iv++)
        if(vmap[iv] == iv && !complete_vertices[iv]) break;
      if(iv == ma->GetNV())//just as a precaution, let us check that it really didnt pitch.
        {
          cout << "Inconsistent data structure. Aborting..." << endl;
          exit(-1);
        }
      for(iv = 0; iv < ma->GetNV(); iv++)
        {
          if(vmap[iv] == iv && !complete_vertices[iv])
            {
              const auto relkt = ktilde[iv] / vrefdt[iv];
              if(relkt < 1e-10) {continue;}
              const auto ttop = tents[latest_tent[iv]]->ttop;
              cout << "v "<<iv<<" tau "<<ttop<<" kt "<<ktilde[iv];
              cout << " rel kt = "<< relkt <<endl;
            }
        }
    }
  delete slabpitcher;
  

  // set lists of internal facets of each element of each tent
  ParallelFor
    (Range(tents),
     [&] (int i)
     {
       Tent & tent = *tents[i];
       TableCreator<int> elfnums_creator(tent.els.Size());

       for ( ; !elfnums_creator.Done(); elfnums_creator++)  {
	 for(int j : Range(tent.els)) {

	   auto fnums = ma->GetElFacets (tent.els[j]);
	   for(int fnum : fnums)
	     if (tent.internal_facets.Pos(fnum) !=
                 tent.internal_facets.ILLEGAL_POSITION)
	       elfnums_creator.Add(j,fnum);
	 }
       }
       tent.elfnums = elfnums_creator.MoveTable();
     });

  // build dependency graph (used by RunParallelDependency)
  TableCreator<int> create_dag(tents.Size());
  for ( ; !create_dag.Done(); create_dag++)
    {
      for (int i : tents.Range())
	for (int d : tents[i]->dependent_tents)
	  create_dag.Add(i, d);
    }
  tent_dependency = create_dag.MoveTable();

  // calculate slope of tents
  ParallelFor
    (Range(tents), [&] (int i)
     {
       LocalHeap slh = lh.Split();
       Tent & tent = *tents[i];

       constexpr auto el_type = EL_TYPE(DIM);
       constexpr int n_vertices = DIM+1; // number of vertices of the current element (simplex)
       ScalarFE<el_type,1> fe; //finite element created for calculating the barycentric coordinates
       IntegrationRule ir(el_type, 0);
       FlatMatrixFixWidth<DIM> dshape_nodal(n_vertices, slh);
       FlatVector<> gradphi_top(DIM, slh), coef_top(n_vertices, slh);
       for (int j : Range(tent.els.Size()))
	 {
	   ElementId ej (VOL, tent.els[j]);
	   const auto vnums = ma->GetElVertices (ej);
	   for (size_t k = 0; k < n_vertices; k++)
	     {
	       auto mapped_vnum = vmap[vnums[k]]; // map periodic vertices
	       auto pos = tent.nbv.Pos(mapped_vnum);
	       if (pos != tent.nbv.ILLEGAL_POSITION)
		 coef_top(k) = tent.nbtime[pos];
	       else
		 coef_top(k) = tent.ttop;
	     }
	   ElementTransformation & trafo = ma->GetTrafo (ej, slh);
	   MappedIntegrationPoint<DIM, DIM> mip(ir[0], trafo);
	   fe.CalcMappedDShape(mip, dshape_nodal);
	   gradphi_top = Trans(dshape_nodal) * coef_top;
	   if(auto norm = L2Norm(gradphi_top); norm > tent.maxslope)
	     tent.maxslope = norm;
	 }
     });
  has_been_pitched = slab_complete;
  return has_been_pitched;
}

template bool TentPitchedSlab::PitchTents<1>(const double, const bool, const double);
template bool TentPitchedSlab::PitchTents<2>(const double, const bool, const double);
template bool TentPitchedSlab::PitchTents<3>(const double, const bool, const double);


double TentPitchedSlab::MaxSlope() const
{
  double maxgrad = 0.0;
  ParallelFor
    (Range(tents),[&] (int i){
      AtomicMax(maxgrad , tents[i]->MaxSlope() );
    });
  return maxgrad;
}


///////////////////// Pitching Algo Routines ///////////////////////////////
TentSlabPitcher::TentSlabPitcher(shared_ptr<MeshAccess> ama, ngstents::PitchingMethod m, Array<int> &avmap) : ma(ama), vertex_refdt(ama->GetNV()), edge_len(ama->GetNEdges()), local_ctau([](const int, const int){return 1.;}), method(m), vmap(avmap) {
  if(method == ngstents::PitchingMethod::EEdgeGrad){ cmax.SetSize(ma->GetNEdges());}
  else {cmax.SetSize(ma->GetNE());}
  cmax = -1;
}


bool TentSlabPitcher::GetReadyVertices(double &adv_factor, bool reset_adv_factor,
                                       const FlatArray<double> &ktilde, const BitArray &complete_vertices,
                                       BitArray &vertex_ready, Array<int> &ready_vertices){

  bool found{false};
  //how many times the adv_factor will be relaxed looking for new vertices
  constexpr int n_attempts = 5;
  vertex_ready = false;
  const double initial_adv_factor = adv_factor;
  for(auto ia = 0; ia < n_attempts; ia++)
    {
      for (auto iv = 0; iv < ma->GetNV(); iv++)
        if(vmap[iv] == iv && !complete_vertices[iv])
          {
            if (ktilde[iv] > adv_factor * vertex_refdt[iv])
              if (!vertex_ready[iv])
                {
                  ready_vertices.Append (iv);
                  vertex_ready.SetBit(iv);
                }
          }
      if(ready_vertices.Size())
        {
          found = true;
          break;
        }
      adv_factor /= 2;
    }
  if(reset_adv_factor)
    adv_factor = initial_adv_factor;
  //the algorithm is most likely stuck
  else if (adv_factor < 0.05) return false;
  return found;
}

void TentSlabPitcher::ComputeVerticesReferenceHeight(const Table<int> &v2v, const Table<int> &v2e, const FlatArray<double> &tau, LocalHeap &lh)
{
  this->vertex_refdt = std::numeric_limits<double>::max();
  for (auto i = 0; i < this->ma->GetNV(); i++)
    if(vmap[i]==i) // non-periodic
      {
        this->vertex_refdt[i] = this->GetPoleHeight(i, tau, v2v[i],v2e[i],lh);
      }
  
}

std::tuple<int,int> TentSlabPitcher::PickNextVertexForPitching(const FlatArray<int> &ready_vertices,
                                                               const FlatArray<double> &ktilde,
                                                               const FlatArray<int> &vertices_level){  
  int minlevel = std::numeric_limits<int>::max();
  int posmin = -1;
  for(auto i = 0; i < ready_vertices.Size(); i++)
    if(vertices_level[ready_vertices[i]] < minlevel)
      {
        minlevel = vertices_level[ready_vertices[i]];
        posmin = i;
      }
  return std::make_tuple(minlevel,posmin);
}

void TentSlabPitcher::UpdateNeighbours(const int vi, const double adv_factor, const Table<int> &v2v,
                                       const Table<int> &v2e, const FlatArray<double> &tau,
                                       const BitArray &complete_vertices, Array<double> &ktilde,
                                       BitArray &vertex_ready, Array<int> &ready_vertices,
                                       LocalHeap &lh){
  for (int nb : v2v[vi])
    {
      nb = vmap[nb]; // map periodic vertices
      if (complete_vertices[nb]) continue;
      const double kt = GetPoleHeight(nb, tau, v2v[nb], v2e[nb],lh);
      ktilde[nb] = kt;
      if (kt > adv_factor * vertex_refdt[nb])
        {
          if (!vertex_ready[nb])
            {
              ready_vertices.Append (nb);
              vertex_ready.SetBit(nb);
            }
        }
      else
        {
          vertex_ready.Clear(nb);
          const auto pos_nb = ready_vertices.Pos(nb);
          if(pos_nb != ready_vertices.ILLEGAL_POSITION)
            {
              ready_vertices.RemoveElement(pos_nb);
            }
        }
    } 
}

template<int DIM>
std::tuple<Table<int>,Table<int>> TentSlabPitcher::InitializeMeshData(LocalHeap &lh, shared_ptr<CoefficientFunction>wavespeed, bool calc_local_ct, const double global_ct)
{
  constexpr auto el_type = EL_TYPE(DIM);//simplex of dimension dim
  constexpr auto n_el_vertices = DIM + 1;//number of vertices of that simplex
  //sets global constant
  this->global_ctau = global_ct;
  BitArray fine_edges(ma->GetNEdges());
  fine_edges.Clear();

  //minimum length of the adjacent edges for each element's vertices
  ArrayMem<double, n_el_vertices> max_edge(n_el_vertices);
  ArrayMem<int,n_el_vertices> v_indices(n_el_vertices);
  //the mesh contains only simplices so only one integration rule is needed
  IntegrationRule ir(el_type, 0);
  for (Ngs_Element el : this->ma->Elements(VOL))
    {
      HeapReset hr(lh);
      max_edge = -1;
      auto ei = ElementId(el);
      ElementTransformation & trafo = this->ma->GetTrafo (ei, lh);
      MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
      const auto wvspd = wavespeed->Evaluate(mip);
      if(method == ngstents::PitchingMethod::EVolGrad)
        {this->cmax[el.Nr()] = wvspd;}


      v_indices = ma->GetElVertices(ei);
      //set all edges belonging to the mesh
      for (int e : el.Edges())
        {
          auto pnts = ma->GetEdgePNums(e);
          auto v1 = pnts[0], v2 = pnts[1];
          if(!fine_edges[e])
            {
              fine_edges.SetBit(e);
              double len = L2Norm (ma-> template GetPoint<DIM>(v1)
                               - ma-> template GetPoint<DIM>(v2));
              edge_len[e] = len;
            }
          if(method == ngstents::PitchingMethod::EEdgeGrad)
            {
              this->cmax[e] = max(this->cmax[e],wvspd);
            }
        }
    }
  //map periodic vertices
  MapPeriodicVertices();
  RemovePeriodicEdges(fine_edges);
  //compute neighbouring data
  TableCreator<int> create_v2e, create_v2v;
  for ( ; !create_v2e.Done(); create_v2e++, create_v2v++)
    {
      for (int e : IntRange (0, ma->GetNEdges()))
        if(fine_edges.Test(e))
          {
            auto vts = ma->GetEdgePNums (e);
            int v1 = vts[0], v2 = vts[1];
            //if v1 (or v2) is not periodic, vmap[v1] == v1
            create_v2v.Add (vmap[v1], v2);
            create_v2e.Add (vmap[v1], e);
            create_v2v.Add (vmap[v2], v1);
            create_v2e.Add (vmap[v2], e);
          }
    }

  TableCreator<int> create_per_verts(ma->GetNV());
  for ( ; !create_per_verts.Done(); create_per_verts++)
    {
      for(auto i : Range(vmap))
        if(vmap[i]!=i)
          create_per_verts.Add(vmap[i],i);
    }

  auto v2v = create_v2v.MoveTable();
  auto v2e = create_v2e.MoveTable();
  per_verts = create_per_verts.MoveTable();
  if(calc_local_ct && DIM > 1)
    {
      local_ctau_table = this->CalcLocalCTau(lh, v2e);
      this->local_ctau = [this](const int v, const int el_or_edge){return local_ctau_table[v][el_or_edge];};
    }
  else
    {
      this->local_ctau = [](const int v, const int el_or_edge){return 1;};
    }
  return std::make_tuple(v2v, v2e);
}

// Get the periodic vertex associated with a primary vertex in periodic case
void TentSlabPitcher::GetVertexElements(int vnr_main, Array<int> & elems) const
{
  ma->GetVertexElements (vnr_main, elems);
  if(per_verts[vnr_main].Size()==0)
    return;
  else
    {
      for(auto per_v : per_verts[vnr_main])
        for(auto elnr : ma->GetVertexElements(per_v))
          elems.Append(elnr);
    }
}

void TentSlabPitcher::GetEdgeElements(int edge, Array<int> & elems) const
{
  ma->GetEdgeElements (edge, elems);

  ArrayMem<int,30> per_edge_els(0);
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_edges = ma->GetPeriodicNodes(NT_EDGE, idnr);
      for (const auto& per_edges : periodic_edges)
        {
          if(per_edges[0] == edge)
            {
              per_edge_els.SetSize(0);
              ma->GetEdgeElements(per_edges[1], per_edge_els);
              for(auto per_el : per_edge_els)
                {elems.Append(per_el);}
            }
        }
    }
}


void TentSlabPitcher::MapPeriodicVertices()
{
  vmap.SetSize(ma->GetNV());
  for (int i : Range(ma->GetNV()))
    vmap[i] = i;
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_nodes = ma->GetPeriodicNodes(NT_VERTEX, idnr);
      for (const auto& per_verts : periodic_nodes)
        vmap[per_verts[1]] = vmap[per_verts[0]];
    }
}


void TentSlabPitcher::RemovePeriodicEdges(BitArray &fine_edges)
{
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_edges = ma->GetPeriodicNodes(NT_EDGE, idnr);
      for (const auto& per_edges : periodic_edges)
        fine_edges.Clear(per_edges[1]);
    }
}



template <int DIM> double VolumeGradientPitcher<DIM>::GetPoleHeight(const int vi, const FlatArray<double> & tau,  FlatArray<int> nbv, FlatArray<int> nbe, LocalHeap & lh) const{
  HeapReset hr(lh);
  constexpr auto el_type = EL_TYPE(DIM);
  //number of vertices of the current element (always the simplex associated to DIM)
  constexpr int n_vertices = DIM+1;
  //finite element created for calculating the barycentric coordinates
  ScalarFE<el_type,1> my_fel;
  // array of all elements containing vertex vi
  ArrayMem<int,30>  els;
  els.SetSize(0);
  this->GetVertexElements(vi, els);

  constexpr double init_pole_height = std::numeric_limits<double>::max();
  double pole_height = init_pole_height;
  //vector containing the advancing front time for each vertex (except vi)
  Vec<n_vertices> coeff_vec(0);
  //gradient of basis functions onthe current element
  FlatMatrixFixWidth<DIM,double> gradphi(n_vertices,lh);
  //indices of el's vertices
  ArrayMem<int,n_vertices> v_indices(n_vertices);
  //numerical tolerance (NOT YET SCALED)
  constexpr double num_tol = std::numeric_limits<double>::epsilon();
  const auto nels = els.Size();
  for (int iel = 0; iel < nels; iel++)
    {
      const auto el = els[iel];
      ElementId ei(VOL,el);
      //c_max^2
      const double c_max_sq = cmax[ei.Nr()] * cmax[ei.Nr()]; 
      //mapping of the current el
      ElementTransformation &trafo = this->ma->GetTrafo(ei, lh);
      //vertices of current el
      v_indices = this->ma->GetElVertices(ei);
      //vi position in current el
      const auto local_vi = v_indices.Pos(vi);
      //integration rule for reference el
      IntegrationRule ir(el_type,1);
      //integration point on deformed element
      MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
      my_fel.CalcMappedDShape(mip,gradphi);

      //sets the coefficient vec
      for (auto k : IntRange(0, v_indices.Size()))
        coeff_vec(k) = tau[v_indices[k]];
      coeff_vec[local_vi] = 0;
      
      /*writing the quadratic eq for tau_v
       \tau_{v}^2 ( \nabla\,\phi_{vi} \cdot \nabla\, \phi_{vi})
                   ^^^^^^^^^^^^^^^^^^alpha^^^^^^^^^^^^^^^^^^^^^^
       +\tau_{v} (\sum_{i \neq v} (\tau_i \nabla,\phi_i \cdot \nabla\,\phi_v))
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^beta^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       + (\sum_{i\neq v}\sum_{j\neq v}\left(\tau_i\tau_j \nabla \phi_i 
``        \cdot \nabla \phi_j\right)-\frac{1}{c}^2)
         ^^^^^^^^^^^^^^^^^^gamma^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

      // square of the norm of gradphi_vi
      const double alpha = InnerProduct(gradphi.Row(local_vi),gradphi.Row(local_vi));
      //since alpha>0 we can scale the equation by alpha
      Vec<DIM> tau_i_grad_i = Trans(gradphi) * coeff_vec;
      const double beta = 2 * InnerProduct(tau_i_grad_i,gradphi.Row(local_vi))/alpha;
      const double gamma = (InnerProduct(tau_i_grad_i, tau_i_grad_i) - 1.0/c_max_sq)/alpha;
      const double delta = beta * beta - 4 * gamma;

      //since sq_delta will always be positive
      //we dont need to consider the solution (-beta-sq_delta)/(2*alpha)
      const double sol = [alpha,beta, gamma, delta,init_pole_height,num_tol](){
        if(delta > num_tol * alpha)//positive delta
          {
            if(beta <= num_tol * alpha)
              return (sqrt(delta)-beta)/2.0;
            else
              return - (2.0 * gamma) / (beta + sqrt(delta));
          }
        else if (delta > -num_tol*alpha)//zero delta
          {
            return -beta/2.0;
          }
        return init_pole_height;//negative delta
      }();
      //the return value is actually the ADVANCE in the current vi
      auto kbar = sol - tau[vi];
      const auto local_ct = local_ctau(vi,iel);
      kbar *= local_ct * global_ctau;
      pole_height = min(pole_height,kbar);
    }
  //check if a real solution to the quadratic equation was found
  //or if the solution is negligible
  //TODO: this checking if a solution was found is quite ugly. how to improve it?
  if( pole_height > 0.75 * init_pole_height) return 0.0;
  else return pole_height * (1 - num_tol);//just to enforce causality
 }

template <int DIM>
Table<double> VolumeGradientPitcher<DIM>::CalcLocalCTau(LocalHeap &lh, const Table<int> &v2e){
  constexpr auto el_type = EL_TYPE(DIM);//simplex of dimension dim
  constexpr auto n_el_vertices = DIM + 1;//number of vertices of that simplex
  
  const auto n_mesh_vertices = ma->GetNV();
  //this table will contain the local mesh-dependent constant
  TableCreator<double> create_local_ctau;
  create_local_ctau.SetSize(n_mesh_vertices);
  create_local_ctau.SetMode(2);
  ArrayMem<int,30>vertex_els(0);
  //just calculating the size of the table
  for(auto vi : IntRange(0, n_mesh_vertices))
    {
      if(vi != vmap[vi]) {continue;}
      this->GetVertexElements(vi,vertex_els);
      const auto n_vert_els = vertex_els.Size();
      for(auto el : IntRange(0, n_vert_els))
        {create_local_ctau.Add(vi,0);}
    }
  create_local_ctau++;// it is in insert mode
  //for a given vertex V in an element E with faces F the constant is calculated as
  //the minimum (over the faces F) ratio between the length of the opposite
  //edge and the biggest edge adjacent to V in F
  //therefore it must be ensured that ctau <=1
  for(auto vi : IntRange(0,n_mesh_vertices))
    {
      if(vi != vmap[vi]){continue;}
      this->GetVertexElements(vi,vertex_els);
      for(auto iel : IntRange(0,vertex_els.Size()))
        {
          HeapReset hr(lh);
          const auto el_num = vertex_els[iel];
          const ElementId ei(VOL,el_num);
          auto faces = ma->GetElFaces(ei);
          double val = 1.0;
          for(auto face : faces)
            {
              auto face_vertices = ma->GetFacePNums(face);
              //check if the face contains the vertex vi
              if(face_vertices.Pos(vi) == face_vertices.ILLEGAL_POSITION)
                {
                  bool found = false;
                  for(auto sl_v : per_verts[vi])
                    {
                      if(face_vertices.Pos(sl_v) == face_vertices.ILLEGAL_POSITION)
                        { found = true;}
                    }
                  if(!found) {continue;}
                }
              
              auto face_edges = ma->GetFaceEdges(face);
              double opposite_edge = -1;
              double max_edge = 0;
              for(auto edge : face_edges)
                {
                  auto pnts = ma->GetEdgePNums(edge);
                  if(vmap[pnts[0]] != vi && vmap[pnts[1]] != vi)
                    {
                      opposite_edge = edge_len[edge];
                    }
                  else
                    {
                      max_edge = max(edge_len[edge],max_edge);
                    }
                }
              val = min(val,opposite_edge/max_edge);
            }
          //val must be <=1
          val = min(val,1.0);
          create_local_ctau.Add(vi,val);
        }
    }

  return create_local_ctau.MoveTable();
}

template <int DIM>
double EdgeGradientPitcher<DIM>::GetPoleHeight(const int vi, const FlatArray<double> & tau, FlatArray<int> nbv, FlatArray<int> nbe, LocalHeap & lh) const{
  double kt = std::numeric_limits<double>::max();
  for (int nb_index : nbv.Range())
    {
      const int nb = vmap[nbv[nb_index]];
      const int edge = nbe[nb_index];
      const double length = edge_len[edge];
      const double c_max = cmax[edge];
      const double local_ct = this->local_ctau(vi,nb_index);
      const double kt1 = tau[nb]-tau[vi]+ global_ctau * local_ct * length/c_max;
      kt = min(kt,kt1);
    }
  constexpr double num_tol = std::numeric_limits<double>::epsilon();
  //ensuring causality (inequality)
  if(kt > num_tol) return kt * (1 - num_tol);
  else return 0.0;
}

template <int DIM>
Table<double> EdgeGradientPitcher<DIM>::CalcLocalCTau(LocalHeap &lh, const Table<int> &v2e)
{
  constexpr auto el_type = EL_TYPE(DIM);//simplex of dimension dim
  constexpr auto n_el_vertices = DIM + 1;//number of vertices of that simplex
  const auto n_mesh_vertices = ma->GetNV();
  //this table will contain the local mesh-dependent constant
  TableCreator<double> create_local_ctau;
  create_local_ctau.SetSize(n_mesh_vertices);
  create_local_ctau.SetMode(2);
  //just calculating the size of the table
  for(auto vi : IntRange(0, n_mesh_vertices))
    {
      if(vi != vmap[vi]) {continue;}
      const auto n_edges_vert = v2e[vi].Size();
      for(auto el : IntRange(0, n_edges_vert))
        {create_local_ctau.Add(vi,0);}
    }
  create_local_ctau++;// it is in insert mode
  
  //used to calculate distance to opposite facet
  ScalarFE<el_type,1> my_fel;
  ArrayMem<int, 30> edge_els(0);
  ArrayMem<int, 30> edge_faces(0);
  //the mesh contains only simplices so only one integration rule is needed
  IntegrationRule ir(el_type, 0);


  //the constant is calculated as the minimum of the projections
  //of the gradient over an edge when the basis functions associated with
  //its vertices are equal to one
  //this constant was developed with the 2D scenario in mind.
  //in 3D, it is thus necessary to scale this projection w.r.t. the
  //projection of the gradient over the respective face
  for(auto vi : IntRange(0, n_mesh_vertices))
    {
      if(vi != vmap[vi]){continue;}
      for(auto edge : v2e[vi])
        {
          //gets the elements that have this edge as a side
          edge_els.SetSize(0);
          
          this->GetEdgeElements(edge, edge_els);
          double val = std::numeric_limits<double>::max();
          //gets the vertices belonging to the edge
          auto pnts = ma->GetEdgePNums(edge);
          auto v1 = pnts[0], v2 = pnts[1];
          //iterate through the elements containing the edge
          for (auto  iel : edge_els)
            {
              HeapReset hr(lh);
              //gradient of basis functions on the current element  
              FlatMatrixFixWidth<DIM,double> gradphi(n_el_vertices,lh);
              const auto ei = ElementId(iel);
              const auto el = ma->GetElement(ei);             

              ElementTransformation &trafo = this->ma->GetTrafo(ei, lh);
              MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
              my_fel.CalcMappedDShape(mip,gradphi);
              //let us test for periodicity support
              auto FindVertexInEl = [el,this](const int v)
              {
                auto &per_vs = this->per_verts[v];
                if(el.Points().Pos(v) != el.Points().ILLEGAL_POSITION)
                  { return el.Points().Pos(v);}
                else
                  {
                    for (const auto &p_v : per_vs)
                      {
                        if(el.Points().Pos(p_v) != el.Points().ILLEGAL_POSITION)
                          { return el.Points().Pos(p_v);}
                      }

                    //we really should not hit this point
                    //it is probably something related to periodicity.
                    throw Exception ("\nngstents error:"
                                     " node numbering inconsistency.\n"
                                     "Please open an issue copying "
                                     "this message.\n");
                    return el.Points().ILLEGAL_POSITION;
                  }
              };
              /*the inner products gradphi.edgevec are identically equal to one so
                there is no need to calculate them*/
              const auto v1_local = FindVertexInEl(vmap[v1]);
              const auto v2_local = FindVertexInEl(vmap[v2]);

              /*
                splitting 2d and 3d code. there is no need to generate that much
                code for 2d
              */
              const auto one_over_max_grad =
                [&]()
                {
                  /*
                    for 2d the projection of the gradient over the (only) face is 
                    always equal to one
                   */
                  if constexpr ( DIM == 2 )
                    {
                      return
                        1.0 / max(L2Norm(gradphi.Row(v1_local)),
                            L2Norm(gradphi.Row(v2_local)));
                    }
                  /*
                   for 3d this is no longer the case
                  */
                  else if constexpr ( DIM == 3 )
                    {
                      edge_faces.SetSize(0);
                      //get all the faces adjacent to the edge
                      ma->GetEdgeFaces(edge, edge_faces);
                      //normal vectors in the REFERENCE element, needed for 3D
                      const auto all_normals = 
                        ElementTopology::GetNormals<DIM>(el_type);
                      Mat<DIM,DIM> inv_jac =  mip.GetJacobianInverse();
                      const double det = fabs(mip.GetJacobiDet());
                      double val = 1;
                      for (auto face : edge_faces)
                        {
                          //local face id
                          const auto face_local = el.Faces().Pos(face);
                          //maybe the current element does not contain this face
                          if(face_local == el.Faces().ILLEGAL_POSITION){continue;}
                          Vec<DIM> normal_ref = all_normals[face_local];
                          //normal vector in the deformed element
                          Vec<DIM> normal = det * Trans(inv_jac) * normal_ref;
                          //the norm of the vector is not unitary
                          const double len_normal = L2Norm(normal);
                          normal /= len_normal;
                          /*
                           this lambda calculates the ratio between
                          the norm of the projection of grad over a face
                          and the norm of grad
                          */
                          auto calc_grad_proj =
                            [&](const int vertex)
                            {
                              Vec<DIM> max_grad_vec = gradphi.Row(vertex);
                              const double norm_grad = L2Norm(gradphi.Row(vertex));
                              max_grad_vec /= norm_grad;
                              const auto tg_grad =
                                L2Norm(Cross(normal, max_grad_vec));
                              return tg_grad / norm_grad;
                            };
                          const double val_face =
                            min(calc_grad_proj(v1_local),
                                calc_grad_proj(v2_local));
                          val = min(val,val_face);
                        }
                      return val;
                    }
                  else//this will never be called by DIM == 1, but anyway
                    {return 1.0;}
                }();
              const auto projGrad = one_over_max_grad /edge_len[edge];
              val = min(val,projGrad);
            }
          create_local_ctau.Add(vi,val);
        }
    }
  return create_local_ctau.MoveTable();
  
}

///////////////////// Output routines //////////////////////////////////////

void TentPitchedSlab::DrawPitchedTentsVTK(string filename)
{
  ofstream out(filename+".vtk");
  Array<Vec<3>> points;
  Array<INT<4>> cells;
  Array<int> level, tentnr;
  int ptcnt = 0;

  for(int i : Range(GetNTents()))
    {
      int firstpt = ptcnt;
      const Tent & tent = GetTent(i);
      Vec<2> pxy = ma->GetPoint<2> (tent.vertex);
      points.Append (Vec<3> (pxy(0), pxy(1), tent.tbot));
      points.Append (Vec<3> (pxy(0), pxy(1), tent.ttop));
      INT<4> tet(ptcnt,ptcnt+1,0,0);
      ptcnt+=2;

      for (int elnr : tent.els)
	{
	  Ngs_Element el = ma->GetElement(ElementId(VOL,elnr));

	  for (int v : el.Vertices())
	    if (v != tent.vertex)
	      {
		pxy = ma->GetPoint<2> (v);
		points.Append (Vec<3> (pxy(0), pxy(1),
                               tent.nbtime[tent.nbv.Pos(v)]));
	      }
	  for (int j = 2; j < 4; j++)
	    tet[j] = ptcnt++;

	  cells.Append(tet);
	}
      for(int j = firstpt; j < ptcnt; j++)
	{
	  level.Append(tent.level);
	  tentnr.Append(i);
	}
    }

  // header
  out << "# vtk DataFile Version 3.0" << endl;
  out << "vtk output" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;


  out << "POINTS " << points.Size() << " float" << endl;
  for (auto p : points)
    out << p << endl;


  out << "CELLS " << cells.Size() << " " << 5 * cells.Size() << endl;
  for (auto c : cells)
    out << 4 <<" " << c << endl;


  out << "CELL_TYPES " << cells.Size() << endl;
  for (auto c : cells)
    out << "10 " << endl;

  out << "CELL_DATA " << cells.Size() << endl;
  out << "POINT_DATA " << points.Size() << endl;

  out << "FIELD FieldData " << 2 << endl;

  out << "tentlevel" << " 1 " << level.Size() << " float" << endl;
  for (auto i : level)
    out << i << " ";
  out << endl;

  out << "tentnumber" << " 1 " << tentnr.Size() << " float" << endl;
  for (auto i : tentnr)
    out << i << " ";
  out << endl;
}


// Used with OpenGL (ngsgui/tents_visualization) and WebGL (webgui).
void TentPitchedSlab::
DrawPitchedTentsGL(Array<int> & tentdata, Array<double> & tenttimes, int & nlevels)
{
  nlevels = 0;
  tentdata.SetAllocSize(4*tents.Size());
  tenttimes.SetAllocSize(4*tents.Size());

  for(int i : Range(GetNTents()))
    {
      const Tent & tent = GetTent(i);
      for(int el : Range(tent.els))
        {
          tentdata.Append(i);
          tentdata.Append(tent.level);
          tentdata.Append(tent.vertex);
          tentdata.Append(tent.els[el]);
          if(tent.level > nlevels)
            nlevels = tent.level;

          if(ma->GetDimension() == 2)
          {
            auto verts = ma->GetElVertices(ElementId(VOL,tent.els[el]));
            for(auto v : verts)
              {
                auto pos = tent.nbv.Pos(vmap[vmap[v]]);
                if (pos != tent.nbv.ILLEGAL_POSITION)
                  tenttimes.Append(tent.nbtime[pos]);
                else
                  tenttimes.Append(tent.tbot);
              }
            tenttimes.Append(tent.ttop);
          }
        }
    }
  nlevels+=1;
}


ostream & operator<< (ostream & ost, const Tent & tent)
{
  ost << "vertex: " << tent.vertex << ", tbot = " << tent.tbot
      << ", ttop = " << tent.ttop << endl;
  ost << "neighbour vertices: " << endl;
  for (size_t k = 0; k < tent.nbv.Size(); k++)
    ost << k << ": " << tent.nbv[k] << " " << tent.nbtime[k] << endl;
  ost << "elements: " << endl << tent.els << endl;
  ost << "internal_facets: " << endl << tent.internal_facets << endl;
  ost << "elfnums: " << endl << tent.elfnums << endl;
  return ost;
}


///////////// TentDataFE ///////////////////////////////////////////////////


TentDataFE::TentDataFE(const Tent & tent, const FESpace & fes, LocalHeap & lh)
  : fei(tent.els.Size(), lh),
    iri(tent.els.Size(), lh),
    miri(tent.els.Size(), lh),
    trafoi(tent.els.Size(), lh),
    mesh_size(tent.els.Size(), lh),
    agradphi_bot(tent.els.Size(), lh),
    agradphi_top(tent.els.Size(), lh),
    adelta(tent.els.Size(), lh),
    felpos(tent.internal_facets.Size(), lh),
    fir(tent.internal_facets.Size(), lh),
    firi(tent.internal_facets.Size(), lh),
    mfiri1(tent.internal_facets.Size(), lh),
    mfiri2(tent.internal_facets.Size(), lh),
    agradphi_botf1(tent.internal_facets.Size(), lh),
    agradphi_topf1(tent.internal_facets.Size(), lh),
    agradphi_botf2(tent.internal_facets.Size(), lh),
    agradphi_topf2(tent.internal_facets.Size(), lh),
    anormals(tent.internal_facets.Size(), lh),
    adelta_facet(tent.internal_facets.Size(), lh)
{
  auto & ma = fes.GetMeshAccess();
  int dim = ma->GetDimension();
  int order = fes.GetOrder();

  size_t ntents = tent.els.Size();
  FlatArray<BaseScalarFiniteElement*> fe_nodal(ntents, lh);
  FlatArray<FlatVector<double>> coef_delta(ntents, lh);
  FlatArray<FlatVector<double>> coef_top(ntents, lh);
  FlatArray<FlatVector<double>> coef_bot(ntents, lh);

  // precompute element data for given tent
  for (size_t i = 0; i < ntents; i++)
    {
      ElementId ei(VOL, tent.els[i]);
      Array<int> dnums;
      fes.GetDofNrs (ei, dnums);
      ranges.Append(IntRange(dnums.Size()) + dofs.Size());
      dofs += dnums;

      fei[i] = &fes.GetFE (ei, lh);
      iri[i] = new (lh) SIMD_IntegrationRule(fei[i]->ElementType(),2*order);
      trafoi[i] = &ma->GetTrafo (ei, lh);
      miri[i] =  &(*trafoi[i]) (*iri[i], lh);

      mesh_size[i] = pow(fabs((*miri[i])[0].GetJacobiDet()[0]),
                         1.0/miri[i]->DimElement());

      auto nipt = miri[i]->Size();
      agradphi_bot[i].AssignMemory(dim, nipt, lh);
      agradphi_top[i].AssignMemory(dim, nipt, lh);
      adelta[i].AssignMemory(nipt, lh);

      switch(dim)
        {
        case 1: fe_nodal[i] = new (lh) ScalarFE<ET_SEGM,1>(); break;
        case 2: fe_nodal[i] = new (lh) ScalarFE<ET_TRIG,1>(); break;
        default: fe_nodal[i] = new (lh) ScalarFE<ET_TET,1>();
        }

      auto ndof = fe_nodal[i]->GetNDof();
      coef_top[i].AssignMemory(ndof, lh);
      coef_bot[i].AssignMemory(ndof, lh);
      auto vnums = ma->GetElVertices(ei);
      auto &vmap = tent.vmap;
      for (size_t k = 0; k < vnums.Size(); k++)
        {
          auto mapped_vnum = vmap[vnums[k]]; // map periodic vertices
          auto pos = tent.nbv.Pos(mapped_vnum);
          if (pos != tent.nbv.ILLEGAL_POSITION)
            coef_bot[i](k) = coef_top[i](k) = tent.nbtime[pos];
          else {
	    coef_bot[i](k) = tent.tbot;
	    coef_top[i](k) = tent.ttop;
	  }
        }
      fe_nodal[i]->EvaluateGrad(*miri[i], coef_top[i], agradphi_top[i]);
      fe_nodal[i]->EvaluateGrad(*miri[i], coef_bot[i], agradphi_bot[i]);

      coef_delta[i].AssignMemory(ndof, lh);
      coef_delta[i] = coef_top[i] - coef_bot[i];
      fe_nodal[i]->Evaluate(*iri[i], coef_delta[i], adelta[i]);
    }
  nd = dofs.Size();

  // precompute facet data for given tent
  for (size_t i = 0; i < tent.internal_facets.Size(); i++)
    {
      INT<2> loc_facetnr;

      ArrayMem<int,2> elnums;
      ArrayMem<int,2> elnums_per;
      ma->GetFacetElements(tent.internal_facets[i], elnums);

      bool periodic_facet = false;
      int facet2;
      if(elnums.Size() < 2)
        {
          facet2 = ma->GetPeriodicFacet(tent.internal_facets[i]);
          if(facet2 != tent.internal_facets[i])
            {
              ma->GetFacetElements (facet2, elnums_per);
              if (elnums_per.Size())
                {
                  periodic_facet = true;
                  elnums.Append(elnums_per[0]);
                }
            }
        }
      
      felpos[i] = INT<2,size_t>(size_t(-1));
      for(int j : Range(elnums.Size()))
        {
          felpos[i][j] = tent.els.Pos(elnums[j]);
          if(felpos[i][j] != size_t(-1))
            {
              auto fnums = ma->GetElFacets (elnums[j]);
              int fnr = tent.internal_facets[i];
              if(periodic_facet)
                {
                  auto pos = fnums.Pos(facet2);
                  if(pos != size_t(-1))
                    fnr = facet2; // change facet nr to periodic
                }
              for (int k : Range(fnums.Size()))
                if (fnums[k] == fnr) loc_facetnr[j] = k;

              auto & trafo = *trafoi[felpos[i][j]];

              auto vnums = ma->GetElVertices (elnums[j]);
              Facet2ElementTrafo transform(trafo.GetElementType(), vnums);

              auto etfacet = ElementTopology::
		GetFacetType (trafo.GetElementType(), loc_facetnr[j]);

              if(j == 0)
                {
		  fir[i] = new (lh) SIMD_IntegrationRule (etfacet, 2*order+1);
		  fir[i]->SetIRX(nullptr); // quick fix to avoid usage of TP elements (slows down)
                }

	      firi[i][j] = &transform(loc_facetnr[j], *fir[i], lh);
	      auto nipt = firi[i][j]->Size();
              if(j == 0)
                {
                  mfiri1[i] = &trafo(*firi[i][j], lh);
                  mfiri1[i]->ComputeNormalsAndMeasure(trafo.GetElementType(),
                                                      loc_facetnr[j]);
                  anormals[i].AssignMemory(dim, nipt, lh);
                  adelta_facet[i].AssignMemory(nipt, lh);
                  agradphi_botf1[i].AssignMemory(dim, nipt, lh);
                  agradphi_topf1[i].AssignMemory(dim, nipt, lh);
		  anormals[i] = Trans(mfiri1[i]->GetNormals());

                  size_t elpos = felpos[i][j];
                  fe_nodal[elpos]->Evaluate(*firi[i][j], coef_delta[elpos],
                                            adelta_facet[i]);
                  fe_nodal[elpos]->EvaluateGrad(*mfiri1[i], coef_bot[elpos],
						agradphi_botf1[i]);
                  fe_nodal[elpos]->EvaluateGrad(*mfiri1[i], coef_top[elpos],
						agradphi_topf1[i]);
                }
              else
                {
                  mfiri2[i] = &trafo(*firi[i][j], lh);
		  agradphi_botf2[i].AssignMemory(dim, nipt, lh);
                  agradphi_topf2[i].AssignMemory(dim, nipt, lh);
                  size_t elpos = felpos[i][j];
                  fe_nodal[elpos]->EvaluateGrad(*mfiri2[i], coef_bot[elpos],
                                                agradphi_botf2[i]);
                  fe_nodal[elpos]->EvaluateGrad(*mfiri2[i], coef_top[elpos],
                                                agradphi_topf2[i]);
                }
            }
        }
    }
}


