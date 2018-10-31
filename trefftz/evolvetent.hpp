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
    template<int D>
    void EvolveTents(int order, shared_ptr<MeshAccess> ma, double wavespeed, double dt, SliceMatrix<> wavefront, double timeshift = 0);

    template<int D>
    Mat<D+1,D+1> TentFaceVerts(Tent* tent, int elnr, shared_ptr<MeshAccess> ma, int top);

    template<int D>
    void TentDmat(Mat<D+1> &Dmat, Mat<D+1> v, int top, double wavespeed);

    template<int D>
    double TentFaceArea( Mat<D+1,D+1> v );

    template<int D>
    Vec<D> TentFaceNormal( Mat<D,D> v, int dir );

    template<int D>
    Vec<D+2,double> TestSolution(Vec<D+1,double> p, double wavespeed);

    template<int D>
    Vector<> EvalBC(const SIMD_MappedIntegrationRule<D,D+1> & mir, double wavespeed, double timeshift);

    template<typename T=double>
    void SwapIfGreater(T& a, T& b);

    template<int D>
    double TentAdiam(Tent* tent, shared_ptr<MeshAccess> ma);
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEvolveTent(py::module m);
#endif // NGS_PYTHON

#endif
