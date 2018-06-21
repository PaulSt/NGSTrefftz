#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>
#include "tents/tents.hpp"

namespace ngcomp
{
  typedef map<netgen::Point3d, netgen::PointIndex> Point2IndexMap;

  void EvolveTents(shared_ptr<MeshAccess> ma);
  shared_ptr<MeshAccess> NgsTPmesh(shared_ptr<MeshAccess> ma, double wavespeed, double dt=1);
  netgen::PointIndex Point2Index(map<netgen::Point3d, netgen::PointIndex> *point2index_map, netgen::Point3d p);
  netgen::PointIndex AddPointUnique(shared_ptr<netgen::Mesh> ngma, Point2IndexMap *pim, netgen::Point3d p);
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m);
#endif // NGS_PYTHON

#endif
