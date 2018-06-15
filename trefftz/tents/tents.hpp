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

  int order;

  int nd;          // total # interior and interface dofs in space
  Array<int> nd_T; // nd_T[k] = # innerdofs of the k-th element of tent
  Array<int> dofs; // all interior and interface dof nums, size(dofs)=nd.
  Array<IntRange> ranges; // ranges[k] IntRange of k-th element of local matrix

  Matrix<> b, mbot, mtop; // nd x nd matrices defined in NumProcTentPitching
  // Matrix<> mstar_inv;
  // Matrix<> propagate;   // propagate * (current u) = u at new time

  Array<Matrix<>> gradphi_bot, gradphi_top;
  Array<AFlatMatrix<>> agradphi_bot, agradphi_top;
  Array<Vector<double>> delta;
  Array<AVector<double>> adelta;
  Array<Vector<>> graddelta;

  Table<Matrix<>> gradphi_facet_bot, gradphi_facet_top;
  Table<Vector<double>> delta_facet;
  Table<AVector<double>> adelta_facet;

  int level;
  int nd_u;                   // num internal dofs
  Array<int> dofs_u;          // internal dofs
  Matrix<> propagate_u;       //
  Array<int> dependent_tents; // these tents depend on me

  class TempTentData *tempdata = nullptr;

  ~Tent ()
  {
    for (auto grad : agradphi_bot)
      free (&grad (0, 0));
    for (auto grad : agradphi_top)
      free (&grad (0, 0));
  }
};

class TempTentData
{
public:
  Array<FiniteElement *> fei;
  Array<SIMD_IntegrationRule *> iri;
  Array<SIMD_BaseMappedIntegrationRule *> miri;
  Array<ElementTransformation *> trafoi;

  TempTentData (int n, LocalHeap &lh)
      : fei (n, lh), iri (n, lh), miri (n, lh), trafoi (n, lh)
  {
    ;
  }
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

  void PitchTents (double dt, double cmax);
  void PitchTents (double dt, double cmax, LocalHeap &lh);

  // Collect tent dofs (inner dofs counted first, then interface dofs),
  // set tent.nd_T, tent.dofs, tent.nd, tent.nd_u, tent.inner_dofs, etc:

  void SetupTents (const shared_ptr<L2HighOrderFESpace> v, LocalHeap &lh);

  int SpacetimeDofs () { return spacetime_dofs; }

  // Export pitched tents into a VTK output file
  void VTKOutputTents (string filename);
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
