#ifndef TENTS_HPP_INCUDED
#define TENTS_HPP_INCUDED

#include <solve.hpp>
#include "paralleldepend.hpp"

using namespace ngsolve;

// A tent is a macroelement consisting of a tentpole erected at a
// central vertex in space and all the space-time tetrahedra attached to
// the tentpole.

// We represent the tent by its projection on space (a vertex patch),
// the central vertex, and the heights (times) of its neighboring
// vertices.

// In this first implementation, we assume that all neighboring
// vertices are either at a current time slice or a new time
// slice. The central vertex is on both time slices.

class Tent
{
public:
  int vertex;           // central vertex
  double tbot, ttop;    // bottom and top times of central vertex
  Array<int> nbv;       // neighbour vertices
  Array<double> nbtime; // time of neighbouring vertices

  Array<int> els;   // all elements in the tent's vertex patch
  Array<int> edges; // all internal facets in the tent's vertex patch
  Table<int>
      elfnums; // elfnums[k] all internal facets of the k-th element of tent

  int nd;          // total # interior and interface dofs in space
  Array<int> nd_T; // nd_T[k] = # innerdofs of the k-th element of tent
  Array<int> dofs; // all interior and interface dof nums, size(dofs)=nd.
  Array<IntRange> ranges; // ranges[k] IntRange of k-th element of local matrix

  Matrix<> b, mbot, mtop; // nd x nd matrices defined in NumProcTentPitching
  // Matrix<> mstar_inv;
  // Matrix<> propagate;   // propagate * (current u) = u at new time

  Array<Matrix<>> gradphi_bot, gradphi_top;
  // Array< AFlatMatrix<> > agradphi_bot, agradphi_top;
  Array<Vector<double>> delta;
  // Array< AVector<double> > adelta;
  Array<Vector<>> graddelta;

  Table<Matrix<>> gradphi_facet_bot, gradphi_facet_top;
  // Table< AFlatMatrix<> > agradphi_facet_bot, agradphi_facet_top;
  Table<Vector<double>> delta_facet;
  // Table< AVector<double> > adelta_facet;

  int level;
  int nd_u;                   // num internal dofs
  Array<int> dofs_u;          // internal dofs
  Matrix<> propagate_u;       //
  Array<int> dependent_tents; // these tents depend on me

  class TempTentData *tempdata = nullptr;

  // global physical times
  double *
      time; // global physical time at vertex, stored in ConservationLaw::gftau
  double timebot; // global physical bottom time at vertex

  void InitTent (shared_ptr<GridFunction> gftau)
  {
    time = &(gftau->GetVector ().FVDouble () (vertex));
    timebot = *time;
  }

  void SetFinalTime () { *time = timebot + (ttop - tbot); }

  ~Tent ()
  {
    // for(auto grad : agradphi_bot)
    //   free(&grad(0,0));
    // for(auto grad : agradphi_top)
    //   free(&grad(0,0));
    // for(auto elgrad : agradphi_facet_bot)
    //   for(auto grad : elgrad)
    //     free(&grad(0,0));
    // for(auto elgrad : agradphi_facet_top)
    //   for(auto grad : elgrad)
    //     free(&grad(0,0));
  }
};

class TempTentData
{
public:
  // element data
  Array<FiniteElement *> fei; // finite elements for all elements in the tent
  Array<SIMD_IntegrationRule *>
      iri; // integration rules for all elements in the tent
  Array<SIMD_BaseMappedIntegrationRule *>
      miri; // mapped integration rules for all elements in the tent
  Array<ElementTransformation *>
      trafoi; // element transformations for all elements in the tent
  Array<double> mesh_size; // mesh size for each element
  Array<FlatMatrix<SIMD<double>>>
      agradphi_bot; // gradient of the old advancing front in the IP's
  Array<FlatMatrix<SIMD<double>>>
      agradphi_top; // gradient of the new advancing front in the IP's
  Array<FlatVector<SIMD<double>>> adelta; // height of the tent in the IP's

  // facet data
  Array<INT<2, size_t>> felpos; // local numbers of the neighbors
  Array<Vec<2, const SIMD_IntegrationRule *>>
      firi; // facet integration rules for all facets in the tent
            // transformed to local coordinated of the
            // neighboring elements
  Array<SIMD_BaseMappedIntegrationRule *>
      mfiri1; // mapped facet integration rules for all facets
  Array<SIMD_BaseMappedIntegrationRule *>
      mfiri2; // mapped facet integration rules for all facets
  Array<FlatMatrix<SIMD<double>>>
      agradphi_botf1; // gradient phi face first and second element
  Array<FlatMatrix<SIMD<double>>> agradphi_topf1;
  Array<FlatMatrix<SIMD<double>>> agradphi_botf2;
  Array<FlatMatrix<SIMD<double>>> agradphi_topf2;

  Array<FlatMatrix<SIMD<double>>> anormals; // normal vectors in the IP's
  // Array<int> facetorder;
  Array<FlatVector<SIMD<double>>>
      adelta_facet; // height of the tent in the IP's

  TempTentData (int n, LocalHeap &lh)
      : fei (n, lh), iri (n, lh), miri (n, lh), trafoi (n, lh)
  {
    ;
  }

  TempTentData (const Tent &tent, const FESpace &fes, const MeshAccess &ma,
                LocalHeap &lh);
};

template <int DIM> class TentPitchedSlab
{

public:
  Array<Tent *> tents; // tents between two time slices
  double dt;           // time step between two time slices
  Table<int> tent_dependency;
  shared_ptr<MeshAccess> ma; // access to base spatial mesh
  int spacetime_dofs;

  TentPitchedSlab (shared_ptr<MeshAccess> ama) : ma (ama) { ; };

  ~TentPitchedSlab ()
  {
    for (auto tent : tents)
      delete tent;
  }

  // Construct tentpitched mesh of slab and tent dependencies using our
  // tent meshing algorithm and initialize all data members of this class:

  double
  GetTentHeight (int vertex, Array<int> &els, FlatArray<int> nbv,
                 Array<double> &tau, Array<double> &cmax, LocalHeap &lh);

  void PitchTents (double dt, double cmax, LocalHeap &lh);
  void
  PitchTents (double dt, shared_ptr<CoefficientFunction> cmax, LocalHeap &lh);
  void PitchTentsGradient (double dt, double cmax, LocalHeap &lh);

  // Collect tent dofs (inner dofs counted first, then interface dofs),
  // set tent.nd_T, tent.dofs, tent.nd, tent.nd_u, tent.inner_dofs, etc:

  void SetupTents (const shared_ptr<L2HighOrderFESpace> v, LocalHeap &lh);

  void CheckTents (const shared_ptr<L2HighOrderFESpace> v, LocalHeap &lh);

  int SpacetimeDofs () { return spacetime_dofs; }

  // Export pitched tents into a VTK output file
  void VTKOutputTents (string filename);

  void GetTentData (Array<int> &tentdata, Array<double> &tentbot,
                    Array<double> &tenttop, int &nlevels);
};

// template class TentPitchedSlab<1>;
// template class TentPitchedSlab<2>;
// template class TentPitchedSlab<3>;

inline ostream &operator<< (ostream &ost, const Tent &tent)
{
  ost << "vertex: " << tent.vertex << ", tbot = " << tent.tbot
      << ", ttop = " << tent.ttop << endl;
  ost << "neighbour vertices: " << endl;
  for (int k = 0; k < tent.nbv.Size (); k++)
    ost << k << ": " << tent.nbv[k] << " " << tent.nbtime[k] << endl;
  ost << "elements: " << endl << tent.els << endl;
  ost << "edges: " << endl << tent.edges << endl;
  ost << "elfnums: " << endl << tent.elfnums << endl;
  ost << "dofs: " << endl << tent.dofs << endl;
  return ost;
}

void VTKOutputTents (shared_ptr<MeshAccess> maptr, Array<Tent *> &tents,
                     string filename);

#endif
