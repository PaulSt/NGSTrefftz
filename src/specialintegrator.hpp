#ifndef SPECIALINTEGRATOR_HPP
#define SPECIALINTEGRATOR_HPP

namespace ngfem
{

  // integrator for \int \lambda(x) \nabla u \nabla v dx
  template <int D>
  class SpaceTimeDG_FFacetBFI : public FacetBilinearFormIntegrator
  {
    double tfix;
    shared_ptr<CoefficientFunction> coef_c;
    shared_ptr<CoefficientFunction> coef_a;
    shared_ptr<CoefficientFunction> coef_sig;
    VorB vb;

  public:
    SpaceTimeDG_FFacetBFI (shared_ptr<CoefficientFunction> acoef_c,
                           shared_ptr<CoefficientFunction> acoef_sig, VorB avb)
        : coef_c (acoef_c), coef_sig (acoef_sig), vb (avb)
    {
      coef_a = make_shared<ConstantCoefficientFunction> (1.0)
               / (acoef_c * acoef_c);
    }

    string Name () const override { return "SpaceTimeDG_FFacetBFI"; }
    VorB VB () const override { return vb; }

    xbool IsSymmetric () const override { return true; }

    void
    CalcElementMatrix (const FiniteElement &fel,
                       const ElementTransformation &eltrans,
                       FlatMatrix<double> elmat, LocalHeap &lh) const override
    {
      throw Exception (
          "SpaceTimeDG_FFacetBFI::CalcElementMatrix - not implemented!");
    }

    void
    CalcFacetMatrix (const FiniteElement &volumefel1, int LocalFacetNr1,
                     const ElementTransformation &eltrans1,
                     FlatArray<int> &ElVertices1,
                     const FiniteElement &volumefel2, int LocalFacetNr2,
                     const ElementTransformation &eltrans2,
                     FlatArray<int> &ElVertices2, FlatMatrix<double> elmat,
                     LocalHeap &lh) const override;

    void
    CalcFacetMatrix (const FiniteElement &volumefel, int LocalFacetNr,
                     const ElementTransformation &eltrans,
                     FlatArray<int> &ElVertices,
                     const ElementTransformation &seltrans,
                     FlatArray<int> &SElVertices, FlatMatrix<double> elmat,
                     LocalHeap &lh) const override;

    void
    ApplyFacetMatrix (const FiniteElement &volumefel1, int LocalFacetNr1,
                      const ElementTransformation &eltrans1,
                      FlatArray<int> &ElVertices1,
                      const FiniteElement &volumefel2, int LocalFacetNr2,
                      const ElementTransformation &eltrans2,
                      FlatArray<int> &ElVertices2, FlatVector<double> elx,
                      FlatVector<double> ely, LocalHeap &lh) const override;

    void
    ApplyFacetMatrix (const FiniteElement &volumefel, int LocalFacetNr,
                      const ElementTransformation &eltrans,
                      FlatArray<int> &ElVertices,
                      const ElementTransformation &seltrans,
                      FlatArray<int> &SElVertices, FlatVector<double> elx,
                      FlatVector<double> ely, LocalHeap &lh) const override;
  };

  // template <int D>
  // class SpaceTimeDG_FFacetLFI : public FacetLinearFormIntegrator
  //{
  // double tfix;
  // shared_ptr<ngcomp::MeshAccess> ma;
  // shared_ptr<CoefficientFunction> coef_c;
  // shared_ptr<CoefficientFunction> coef_a;
  // shared_ptr<CoefficientFunction> coef_sig;
  // VorB vb;

  // public:
  // SpaceTimeDG_FFacetLFI (shared_ptr<ngcomp::MeshAccess> ama,
  // shared_ptr<CoefficientFunction> acoef_c,
  // shared_ptr<CoefficientFunction> acoef_sig,
  // VorB avb)
  //: ma (ama), coef_c (acoef_c), coef_sig (acoef_sig), vb (avb)
  //{
  // coef_a
  //= make_shared<ConstantCoefficientFunction> (1) / (acoef_c * acoef_c);
  //}

  // string Name () const override { return "SpaceTimeDG_FFacetLFI"; }
  // VorB VB () const override { return vb; }

  ////bool BoundaryForm () const override { return false; }

  //// Calculates the right hand side element vector
  // void
  // CalcElementVector (const FiniteElement &fel,
  // const ElementTransformation &eltrans,
  // FlatVector<double> elvec, LocalHeap &lh) const override;
  //};
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSpecialIntegrator (py::module m);
#endif

#endif
