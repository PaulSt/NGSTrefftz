#pragma once

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
// #include <variant>

#ifndef FILE_INTEGRATORCFHPP
#include <integratorcf.hpp> //
#define FILE_INTEGRATORCFHPP
#endif

// #include <variant> //
// #include <memory>  //

#include <comp.hpp>
using namespace ngcomp;

namespace ngfem
{

  class BoxDifferentialSymbol;
  class BoxIntegral : public Integral
  {
    double reference_box_length = 0.5;

  public:
    BoxIntegral (shared_ptr<CoefficientFunction> _cf, DifferentialSymbol _dx,
                 double _reference_box_length);
    BoxIntegral (shared_ptr<CoefficientFunction> _cf,
                 shared_ptr<BoxDifferentialSymbol> _dx);
    virtual ~BoxIntegral () {}

    template <typename TSCAL, int D>
    TSCAL T_BoxIntegrate (const ngcomp::MeshAccess &ma,
                          FlatVector<TSCAL> element_wise);

    virtual double Integrate (const ngcomp::MeshAccess &ma,
                              FlatVector<double> element_wise) override;

    virtual Complex Integrate (const ngcomp::MeshAccess &ma,
                               FlatVector<Complex> element_wise) override;

    virtual shared_ptr<BilinearFormIntegrator>
    MakeBilinearFormIntegrator () const override;
    virtual shared_ptr<LinearFormIntegrator>
    MakeLinearFormIntegrator () const override;

    virtual shared_ptr<Integral>
    CreateSameIntegralType (shared_ptr<CoefficientFunction> _cf) override
    {
      return make_shared<BoxIntegral> (_cf, dx, reference_box_length);
    }
  };

  class BoxDifferentialSymbol : public DifferentialSymbol
  {
    double reference_box_length = 0.5;
    friend class BoxIntegral;

  public:
    BoxDifferentialSymbol (const BoxDifferentialSymbol &bds)
        : DifferentialSymbol (bds),
          reference_box_length (bds.reference_box_length)
    {
      ;
    }

    BoxDifferentialSymbol (double _reference_box_length = 0.5)
        : DifferentialSymbol (VOL),
          reference_box_length (_reference_box_length)
    {
      ;
    }

    virtual ~BoxDifferentialSymbol () { ; }

    virtual shared_ptr<Integral>
    MakeIntegral (shared_ptr<CoefficientFunction> cf) const
    {
      return make_shared<BoxIntegral> (cf, *this, reference_box_length);
    }
  };

}

class BoxLinearFormIntegrator : public SymbolicLinearFormIntegrator
{
public:
  double reference_box_length = 0.5;
  BoxLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB vb = VOL,
                           double _reference_box_length = 0.5);

  virtual string Name () const override { return string ("BoxInt-LFI"); }

  virtual void
  CalcElementVector (const FiniteElement &fel,
                     const ElementTransformation &trafo,
                     FlatVector<double> elvec, LocalHeap &lh) const override;

  virtual void
  CalcElementVector (const FiniteElement &fel,
                     const ElementTransformation &trafo,
                     FlatVector<Complex> elvec, LocalHeap &lh) const override;

  template <int D, typename SCAL>
  void T_CalcElementVector (const FiniteElement &fel,
                            const ElementTransformation &trafo,
                            FlatVector<SCAL> elvec, LocalHeap &lh) const;
};

class BoxBilinearFormIntegrator : public SymbolicBilinearFormIntegrator
{
  double reference_box_length = 0.5;

public:
  BoxBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf,
                             VorB vb = VOL,
                             double _reference_box_length = 0.5);

  virtual xbool IsSymmetric () const override { return maybe; }
  virtual string Name () const override { return string ("BoxInt-BFI"); }

  virtual void
  CalcElementMatrix (const FiniteElement &fel,
                     const ElementTransformation &trafo,
                     FlatMatrix<double> elmat, LocalHeap &lh) const override;

  virtual void
  CalcElementMatrixAdd (const FiniteElement &fel,
                        const ElementTransformation &trafo,
                        FlatMatrix<double> elmat, bool &symmetric_so_far,
                        LocalHeap &lh) const override;

  virtual void
  CalcElementMatrixAdd (const FiniteElement &fel,
                        const ElementTransformation &trafo,
                        FlatMatrix<Complex> elmat, bool &symmetric_so_far,
                        LocalHeap &lh) const override;

  template <int D, typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void
  T_CalcElementMatrixAdd (const FiniteElement &fel,
                          const ElementTransformation &trafo,
                          FlatMatrix<SCAL_RES> elmat, LocalHeap &lh) const;

  template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void
  T_CalcElementMatrixEBAdd (const FiniteElement &fel,
                            const ElementTransformation &trafo,
                            FlatMatrix<SCAL_RES> elmat, LocalHeap &lh) const;

  virtual void CalcLinearizedElementMatrix (const FiniteElement &fel,
                                            const ElementTransformation &trafo,
                                            FlatVector<double> elveclin,
                                            FlatMatrix<double> elmat,
                                            LocalHeap &lh) const override;

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void
  T_CalcLinearizedElementMatrixEB (const FiniteElement &fel,
                                   const ElementTransformation &trafo,
                                   FlatVector<SCAL> elveclin,
                                   FlatMatrix<SCAL> elmat, LocalHeap &lh) const
  {
    throw Exception ("BoxBilinearFormIntegrator::T_"
                     "CalcLinearizedElementMatrixEB not yet implemented");
  }

  virtual void
  ApplyElementMatrix (const FiniteElement &fel,
                      const ElementTransformation &trafo,
                      const FlatVector<double> elx, FlatVector<double> ely,
                      void *precomputed, LocalHeap &lh) const override;

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void
  T_ApplyElementMatrixEB (const FiniteElement &fel,
                          const ElementTransformation &trafo,
                          const FlatVector<SCAL> elx, FlatVector<SCAL> ely,
                          void *precomputed, LocalHeap &lh) const
  {
    throw Exception ("BoxBilinearFormIntegrator::T_ApplyElementMatrixEB not "
                     "yet implemented");
  }
};

#include <python_comp.hpp>

#ifdef NGS_PYTHON
void ExportBoxIntegral (py::module m);
#endif
