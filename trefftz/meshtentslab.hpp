#ifndef FILE_MESHTENTSLAB_HPP
#define FILE_MESHTENTSLAB_HPP

#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>
#include "tents/tents.hpp"

namespace ngcomp
{
  typedef map<netgen::Point3d, netgen::PointIndex> Point2IndexMap;

  shared_ptr<MeshAccess>
  NgsTPmesh (shared_ptr<MeshAccess> ma, double wavespeed, double dt);
  netgen::PointIndex
  Point2Index (map<netgen::Point3d, netgen::PointIndex> *point2index_map,
               netgen::Point3d p);
  netgen::PointIndex AddPointUnique (shared_ptr<netgen::Mesh> ngma,
                                     Point2IndexMap *pim, netgen::Point3d p);
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMeshTentSlab (py::module m);
#endif // NGS_PYTHON

#endif
