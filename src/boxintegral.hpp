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
  enum BOXTYPE
  {
    DEFAULT = -1,
    BOX = 0,
    BALL = 1
  };

  class BoxDifferentialSymbol;
  class BoxIntegral : public Integral
  {
  public:
    double box_length = 0.5;
    bool scale_with_elsize = false;
    BOXTYPE boxtype = DEFAULT;

    BoxIntegral (shared_ptr<CoefficientFunction> _cf, DifferentialSymbol _dx,
                 double _box_length, bool _scale_with_elsize,
                 BOXTYPE _boxtype);
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
      return make_shared<BoxIntegral> (_cf, dx, box_length, scale_with_elsize,
                                       boxtype);
    }
  };

  class BoxDifferentialSymbol : public DifferentialSymbol
  {
  public:
    friend class BoxIntegral;
    double box_length = 0.5;
    bool scale_with_elsize = false;
    BOXTYPE boxtype = DEFAULT;

    BoxDifferentialSymbol (const BoxDifferentialSymbol &bds)
        : DifferentialSymbol (bds), box_length (bds.box_length),
          scale_with_elsize (bds.scale_with_elsize), boxtype (bds.boxtype)
    {
      ;
    }

    BoxDifferentialSymbol (double _box_length = 0.5,
                           bool _scale_with_elsize = false,
                           BOXTYPE _boxtype = DEFAULT)
        : DifferentialSymbol (VOL), box_length (_box_length),
          scale_with_elsize (_scale_with_elsize), boxtype (_boxtype)
    {
      ;
    }

    virtual ~BoxDifferentialSymbol () { ; }

    virtual shared_ptr<Integral>
    MakeIntegral (shared_ptr<CoefficientFunction> cf) const
    {
      return make_shared<BoxIntegral> (cf, *this, box_length,
                                       scale_with_elsize, boxtype);
    }
  };

}

class BoxLinearFormIntegrator : public SymbolicLinearFormIntegrator
{
  double box_length = 0.5;
  bool scale_with_elsize = false;
  BOXTYPE boxtype = DEFAULT;

public:
  BoxLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB vb = VOL,
                           double _box_length = 0.5,
                           bool _scale_with_elsize = false,
                           BOXTYPE _boxtype = DEFAULT);

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
  double box_length = 0.5;
  bool scale_with_elsize = false;
  BOXTYPE boxtype = DEFAULT;

public:
  BoxBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf,
                             VorB vb = VOL, double _box_length = 0.5,
                             bool _scale_with_elsize = false,
                             BOXTYPE _boxtype = DEFAULT);

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
  void T_CalcLinearizedElementMatrixEB (const FiniteElement &,
                                        const ElementTransformation &,
                                        FlatVector<SCAL>, FlatMatrix<SCAL>,
                                        LocalHeap &) const
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
  T_ApplyElementMatrixEB (const FiniteElement &, const ElementTransformation &,
                          const FlatVector<SCAL>, FlatVector<SCAL>, void *,
                          LocalHeap &) const
  {
    throw Exception ("BoxBilinearFormIntegrator::T_ApplyElementMatrixEB not "
                     "yet implemented");
  }
};

#include <python_comp.hpp>

#ifdef NGS_PYTHON
void ExportBoxIntegral (py::module m);
#endif
