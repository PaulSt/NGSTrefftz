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

    template<int D>
    void EvolveTents(shared_ptr<MeshAccess> ma, double wavespeed, double dt, Vector<double> wavefront);

    template<int D>
    double AreaTentFace(Tent* tent, int elnr, shared_ptr<MeshAccess> ma);

    template<int D>
    Vec<D+1> VertexTimesTentFace(Tent* tent, const INT<D+1>& verts);
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m);
#endif // NGS_PYTHON

#endif
