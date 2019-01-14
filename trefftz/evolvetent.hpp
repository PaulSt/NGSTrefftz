#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <fem.hpp>
#include <multigrid.hpp>
#include "tents/tents.hpp"
#include "trefftzwavefe.hpp"

namespace ngcomp
{
    template<int D>
    void EvolveTents(int order, shared_ptr<MeshAccess> ma, double wavespeed, double dt, SliceMatrix<> wavefront, 
            double timeshift,  shared_ptr<CoefficientFunction> bddatum);

    template<int D>
    void CalcTentEl(int elnr, Tent* tent, TrefftzWaveFE<D+1> tel, shared_ptr<MeshAccess> ma, SliceMatrix<double> &wavefront, 
            SIMD_IntegrationRule &sir, LocalHeap &slh, FlatMatrix<> &elmat, FlatVector<> &elvec);

    template<int D>
    void CalcTentBndEl(int surfel, Tent* tent, TrefftzWaveFE<D+1> tel, shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> bddatum, 
                       double timeshift, SIMD_IntegrationRule &sir, LocalHeap &slh, FlatMatrix<> &elmat, FlatVector<> &elvec);

    template<int D>
    void CalcTentElEval(int elnr, Tent* tent, TrefftzWaveFE<D+1> tel, shared_ptr<MeshAccess> ma , SliceMatrix<> &wavefront, SIMD_IntegrationRule &sir, LocalHeap &slh, FlatVector<> &sol);

    template<int D>
    Mat<D+1,D+1> TentFaceVerts(Tent* tent, int elnr, shared_ptr<MeshAccess> ma, int top);

    template<int D>
    void TentDmat(Mat<D+1> &Dmat, Mat<D+1> v, int top, double wavespeed);

    template<int D>
    double TentFaceArea( Mat<D+1,D+1> v );

    template<int D>
    Vec<D> TentFaceNormal( Mat<D,D> v, int dir );

    template<int D>
    Matrix<> MakeWavefront(int order, shared_ptr<MeshAccess> ma, double time, shared_ptr<CoefficientFunction> bddatum);

    template<int D>
    double Error(int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront, Matrix<> wavefront_corr);

    template<int D>
    double Energy(int order, shared_ptr<MeshAccess> ma, Matrix<> wavefront, double wavenumber);

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
