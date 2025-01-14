#ifndef SPECIALINTEGRATOR_HPP
#define SPECIALINTEGRATOR_HPP

#include <symbolicintegrator.hpp>

#include "ngsttd.hpp"

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
    CalcElementMatrix (const FiniteElement &, const ElementTransformation &,
                       FlatMatrix<double>, LocalHeap &) const override
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

  template <int D>
  class SpaceTimeDG_FFacetLFI : public FacetLinearFormIntegrator
  {
    double tfix;
    shared_ptr<ngcomp::MeshAccess> ma;
    shared_ptr<CoefficientFunction> gfuh;
    shared_ptr<CoefficientFunction> gfduh;
    shared_ptr<CoefficientFunction> coef_c;
    shared_ptr<CoefficientFunction> coef_a;
    shared_ptr<CoefficientFunction> coef_sig;
    VorB vb;

  public:
    SpaceTimeDG_FFacetLFI (shared_ptr<ngcomp::MeshAccess> ama,
                           shared_ptr<CoefficientFunction> agfuh,
                           shared_ptr<CoefficientFunction> agfduh,
                           shared_ptr<CoefficientFunction> acoef_c,
                           shared_ptr<CoefficientFunction> acoef_sig, VorB avb)
        : ma (ama), gfuh (agfuh), gfduh (agfduh), coef_c (acoef_c),
          coef_sig (acoef_sig), vb (avb)
    {
      coef_a
          = make_shared<ConstantCoefficientFunction> (1) / (acoef_c * acoef_c);
    }

    string Name () const override { return "SpaceTimeDG_FFacetLFI"; }
    VorB VB () const override { return vb; }

    // bool BoundaryForm () const override { return false; }
    void
    CalcFacetVector (const FiniteElement &volumefel1, int LocalFacetNr1,
                     const ElementTransformation &eltrans1,
                     FlatArray<int> &ElVertices1,
                     const FiniteElement &volumefel2, int LocalFacetNr2,
                     const ElementTransformation &eltrans2,
                     FlatArray<int> &ElVertices2, FlatVector<double> elvec,
                     LocalHeap &lh) const override;

    void
    CalcFacetVector (const FiniteElement &volumefel, int LocalFacetNr,
                     const ElementTransformation &eltrans,
                     FlatArray<int> &ElVertices,
                     const ElementTransformation &seltrans,
                     FlatVector<double> elvec, LocalHeap &lh) const override;
  };

  class SymbolicFFacetBilinearFormIntegrator
      : public FacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction *> trial_proxies, test_proxies;
    Array<CoefficientFunction *> gridfunction_cfs;
    Array<CoefficientFunction *> cache_cfs;
    Array<int> trial_cum, test_cum; // cumulated dimension of proxies
    VorB vb;
    bool element_boundary;
    bool neighbor_testfunction;
    Array<shared_ptr<CoefficientFunction>>
        dcf_dtest; // derivatives by test-functions
  public:
    NGST_DLL
    SymbolicFFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf,
                                          VorB avb, bool aelement_boundary);

    virtual VorB VB () const { return vb; }
    virtual bool BoundaryForm () const { return vb == BND; }
    virtual xbool IsSymmetric () const { return maybe; }

    virtual DGFormulation GetDGFormulation () const
    {
      return DGFormulation (neighbor_testfunction, element_boundary);
    }

    NGST_DLL virtual void
    CalcFacetMatrix (const FiniteElement &volumefel1, int LocalFacetNr1,
                     const ElementTransformation &eltrans1,
                     FlatArray<int> &ElVertices1,
                     const FiniteElement &volumefel2, int LocalFacetNr2,
                     const ElementTransformation &eltrans2,
                     FlatArray<int> &ElVertices2, FlatMatrix<double> elmat,
                     LocalHeap &lh) const;

    NGST_DLL virtual void
    CalcFacetMatrix (const FiniteElement &volumefel, int LocalFacetNr,
                     const ElementTransformation &eltrans,
                     FlatArray<int> &ElVertices,
                     const ElementTransformation &seltrans,
                     FlatArray<int> &SElVertices, FlatMatrix<double> elmat,
                     LocalHeap &lh) const;

    NGST_DLL virtual void
    CalcFacetMatrix (const FiniteElement &volumefel1, int LocalFacetNr1,
                     const ElementTransformation &eltrans1,
                     FlatArray<int> &ElVertices1,
                     const FiniteElement &volumefel2, int LocalFacetNr2,
                     const ElementTransformation &eltrans2,
                     FlatArray<int> &ElVertices2, FlatMatrix<Complex> elmat,
                     LocalHeap &lh) const;

    NGST_DLL virtual void
    CalcFacetMatrix (const FiniteElement &volumefel, int LocalFacetNr,
                     const ElementTransformation &eltrans,
                     FlatArray<int> &ElVertices,
                     const ElementTransformation &seltrans,
                     FlatArray<int> &SElVertices, FlatMatrix<Complex> elmat,
                     LocalHeap &lh) const;

    NGST_DLL virtual void
    CalcLinearizedFacetMatrix (const FiniteElement &, int,
                               const ElementTransformation &, FlatArray<int> &,
                               const ElementTransformation &, FlatArray<int> &,
                               FlatVector<double>, FlatMatrix<double>,
                               LocalHeap &) const
    {
      throw Exception (
          "SymbolicFFacetBFI::CalcLinearizedFacetMatrix not implemented");
    }

    NGST_DLL virtual void
    ApplyFacetMatrix (const FiniteElement &, int,
                      const ElementTransformation &, FlatArray<int> &,
                      const FiniteElement &, int,
                      const ElementTransformation &, FlatArray<int> &,
                      FlatVector<double>, FlatVector<double>,
                      LocalHeap &) const
    {
      throw Exception ("SymbolicFFacetBFI::ApplyFacetMatrix not implemented");
    }

    NGST_DLL virtual void
    CalcTraceValues (const FiniteElement &, int, const ElementTransformation &,
                     FlatArray<int> &, FlatVector<double> &,
                     FlatVector<double>, LocalHeap &) const
    {
      throw Exception ("SymbolicFFacetBFI::CalcTraceValues not implemented");
    }

    NGST_DLL virtual void
    ApplyFromTraceValues (const FiniteElement &, int,
                          const ElementTransformation &, FlatArray<int> &,
                          FlatVector<double>, FlatVector<double>,
                          FlatVector<double>, LocalHeap &) const
    {
      throw Exception (
          "SymbolicFFacetBFI::ApplyFromTraceValues not implemented");
    }

    NGST_DLL virtual void
    ApplyFacetMatrix (const FiniteElement &, int,
                      const ElementTransformation &, FlatArray<int> &,
                      const ElementTransformation &, FlatArray<int> &,
                      FlatVector<double>, FlatVector<double>,
                      LocalHeap &) const
    {
      throw Exception ("SymbolicFFacetBFI::ApplyFacetMatrix not implemented");
    }

  private:
    template <typename TSCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix (const FiniteElement &volumefel1, int LocalFacetNr1,
                            const ElementTransformation &eltrans1,
                            FlatArray<int> &ElVertices1,
                            const FiniteElement &volumefel2, int LocalFacetNr2,
                            const ElementTransformation &eltrans2,
                            FlatArray<int> &ElVertices2,
                            FlatMatrix<TSCAL> elmat, LocalHeap &lh) const;

    template <typename TSCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix (const FiniteElement &volumefel, int LocalFacetNr,
                            const ElementTransformation &eltrans,
                            FlatArray<int> &ElVertices,
                            const ElementTransformation &seltrans,
                            FlatArray<int> &SElVertices,
                            FlatMatrix<TSCAL> elmat, LocalHeap &lh) const;
  };

  class SymbolicFFacetLinearFormIntegrator : public FacetLinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction *> proxies;
    Array<CoefficientFunction *> cache_cfs;
    Array<int> test_cum; // cumulated dimension of proxies
    VorB vb;             // only BND supported by now
    // bool element_boundary;  /// not needed (by now ???)
    IntegrationRule ir;           // if non-empty use this integration-rule
    SIMD_IntegrationRule simd_ir; // if non-empty use this integration-rule

  public:
    SymbolicFFacetLinearFormIntegrator (shared_ptr<CoefficientFunction> acf,
                                        VorB avb);

    virtual VorB VB () const override { return vb; }
    virtual bool BoundaryForm () const override { return vb == BND; }

    NGST_DLL void
    CalcFacetVector (const FiniteElement &volumefel1, int LocalFacetNr1,
                     const ElementTransformation &eltrans1,
                     FlatArray<int> &ElVertices1,
                     const FiniteElement &volumefel2, int LocalFacetNr2,
                     const ElementTransformation &eltrans2,
                     FlatArray<int> &ElVertices2, FlatVector<double> elvec,
                     LocalHeap &lh) const override;

    NGST_DLL void
    CalcFacetVector (const FiniteElement &volumefel1, int LocalFacetNr1,
                     const ElementTransformation &eltrans1,
                     FlatArray<int> &ElVertices1,
                     const FiniteElement &volumefel2, int LocalFacetNr2,
                     const ElementTransformation &eltrans2,
                     FlatArray<int> &ElVertices2, FlatVector<Complex> elvec,
                     LocalHeap &lh) const override;

    NGST_DLL void
    CalcFacetVector (const FiniteElement &volumefel, int LocalFacetNr,
                     const ElementTransformation &eltrans,
                     FlatArray<int> &ElVertices,
                     const ElementTransformation &seltrans,
                     FlatVector<double> elvec, LocalHeap &lh) const override;

    NGST_DLL void
    CalcFacetVector (const FiniteElement &volumefel, int LocalFacetNr,
                     const ElementTransformation &eltrans,
                     FlatArray<int> &ElVertices,
                     const ElementTransformation &seltrans,
                     FlatVector<Complex> elvec, LocalHeap &lh) const override;

  private:
    template <typename TSCAL>
    void
    T_CalcFacetVector (const FiniteElement &fel1, int LocalFacetNr1,
                       const ElementTransformation &trafo1,
                       FlatArray<int> &ElVertices1, const FiniteElement &fel2,
                       int LocalFacetNr2, const ElementTransformation &trafo2,
                       FlatArray<int> &ElVertices2, FlatVector<TSCAL> elvec,
                       LocalHeap &lh) const;

    template <typename TSCAL>
    void T_CalcFacetVector (const FiniteElement &volumefel, int LocalFacetNr,
                            const ElementTransformation &eltrans,
                            FlatArray<int> &ElVertices,
                            const ElementTransformation &seltrans,
                            FlatVector<TSCAL> elvec, LocalHeap &lh) const;
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSpecialIntegrator (py::module m);
#endif

#endif
