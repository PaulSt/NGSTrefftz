#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP

#include "scalarmappedfe.hpp"

#include <fespace.hpp>

/// Denotes the types of equations supported by the TrefftzFESpace.
enum class EqType
{
  fowave,         /// for the first order acoustic wave equation
  foqtwave,       /// for the quasi-Trefftz space related to fowave
  wave,           /// for the second order acoustic wave equation
  qtwave,         /// for the quasi-Trefftz space related to wave
  fowave_reduced, /// derivatives of Trefftz second order wave equation basis
                  /// without constants
  heat,           /// caloric polynomials
  qtheat,         /// quasi-Trefftz space for the heat equation
  laplace,        /// for Laplace equation
  qtelliptic,     /// quasi-Trefftz space for a general elliptic problem
  helmholtz,      /// plane waves for the Helmholtz equation
  helmholtzconj,  /// complex conjugate of Helmholtz
};

namespace ngcomp
{
  using namespace std;
  using namespace ngfem;

  class PolBasis
  {
  protected:
    int order;

  public:
    PolBasis () { ; }
    PolBasis (int aorder) : order (aorder) { ; }
    virtual ~PolBasis () { ; }

    virtual void SetRHS (shared_ptr<CoefficientFunction>)
    {
      throw Exception ("SetRHS not implemented for this basis");
    }
    virtual void
    GetParticularSolution (Vec<1>, Vec<1>, FlatVector<>, LocalHeap &)
    {
      throw Exception ("GetParticularSolution not implemented for this basis");
    }
    virtual void
    GetParticularSolution (Vec<2>, Vec<2>, FlatVector<>, LocalHeap &)
    {
      throw Exception ("GetParticularSolution not implemented for this basis");
    }
    virtual void
    GetParticularSolution (Vec<3>, Vec<3>, FlatVector<>, LocalHeap &)
    {
      throw Exception ("GetParticularSolution not implemented for this basis");
    }
    template <int D>
    void ComputeDerivs (int order, shared_ptr<CoefficientFunction> acoeff,
                        Vector<shared_ptr<CoefficientFunction>> &ders)
    {
      const int ndiffs = (BinCoeff (D + order, order));
      ders.SetSize (ndiffs);

      shared_ptr<CoefficientFunction> coeff = acoeff;
      shared_ptr<CoefficientFunction> coeffx = acoeff;
      shared_ptr<CoefficientFunction> coeffxy = acoeff;

      switch (D)
        {
        case 0:
          break;
        case 1:
          for (int i = 0; i <= order; i++)
            {
              ders (i) = coeffx;
              coeffx = coeffx->Diff (
                  MakeCoordinateCoefficientFunction (0).get (),
                  make_shared<ConstantCoefficientFunction> (1));
            }
          break;
        case 2:
          for (int i = 0; i <= order; i++)
            {
              for (int j = 0; j <= order - i; j++)
                {
                  int iix = IndexMap2<2> (Vec<2, int>{ j, i }, order);
                  ders (iix) = coeffx;
                  coeffx = coeffx->Diff (
                      MakeCoordinateCoefficientFunction (0).get (),
                      make_shared<ConstantCoefficientFunction> (1));
                }
              coeff
                  = coeff->Diff (MakeCoordinateCoefficientFunction (1).get (),
                                 make_shared<ConstantCoefficientFunction> (1));
              coeffx = coeff;
            }
          break;
        case 3:
          for (int i = 0; i <= order; i++)
            {
              for (int j = 0; j <= order - i; j++)
                {
                  for (int k = 0; k <= order - i - j; k++)
                    {
                      int iix = IndexMap2<3> (Vec<3, int>{ k, j, i }, order);
                      ders (iix) = coeffxy;
                      coeffxy = coeffxy->Diff (
                          MakeCoordinateCoefficientFunction (0).get (),
                          make_shared<ConstantCoefficientFunction> (1));
                    }
                  coeffx = coeffx->Diff (
                      MakeCoordinateCoefficientFunction (1).get (),
                      make_shared<ConstantCoefficientFunction> (1));
                  coeffxy = coeffx;
                }
              coeff
                  = coeff->Diff (MakeCoordinateCoefficientFunction (2).get (),
                                 make_shared<ConstantCoefficientFunction> (1));
              coeffx = coeff;
              coeffxy = coeff;
            }
          break;
        }
    }

    template <int D> static int IndexMap2 (Vec<D, int> index, int ord)
    {
      int sum = 0;
      int temp_size = 0;
      for (int d = 0; d < D; d++)
        {
          for (int p = 0; p < index (d); p++)
            {
              sum += BinCoeff (D - 1 - d + ord - p - temp_size,
                               ord - p - temp_size);
            }
          temp_size += index (d);
        }
      return sum;
    }
    // TODO: virtual static CSR Basis(int ord, int basistype = 0, int fowave =
    // 0);
  };

  class TrefftzFESpace : public FESpace
  {
    int D;
    int nel;
    int local_ndof;
    double coeff_const = 1;
    EqType eqtype = EqType::wave;
    int useshift = 1;
    int usescale = 1;
    int basistype = 0;
    shared_ptr<CoefficientFunction> coeffA = nullptr;
    shared_ptr<CoefficientFunction> coeffB = nullptr;
    shared_ptr<CoefficientFunction> coeffC = nullptr;

    CSR basismat;
    Vector<CSR> basismats;
    PolBasis *basis = nullptr;

    // The following functions are helper functions, that capture behavior
    // which strongly depends on the dimension and the eqtype of the
    // TrefftzFESpace.

    /// @returns the local number of dofs for this TrefftzFESpace
    int calcLocalNdofs () const;
    template <int Dim> void setupEvaluators ();
    /// templated version of UpdateBasis
    template <int Dim> void basisUpdate ();
    /// templated version of GetFe
    template <int Dim>
    FiniteElement &TGetFE (ElementId ei, Allocator &alloc) const;

  public:
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
                    bool checkflags = false);
    ~TrefftzFESpace ();
    void Update () override;
    void UpdateCouplingDofArray () override;
    void SetCoeff (double acoeff_const);
    void SetCoeff (shared_ptr<CoefficientFunction> acoeffA,
                   shared_ptr<CoefficientFunction> acoeffB = nullptr,
                   shared_ptr<CoefficientFunction> acoeffC = nullptr);

    shared_ptr<GridFunction>
    GetParticularSolution (shared_ptr<CoefficientFunction> acoeffF);
    string GetClassName () const override { return "trefftz"; }
    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;
    static DocInfo GetDocu ();

  protected:
    void UpdateBasis ();
    template <int D>
    double ElSize (ElementId ei, Vec<D> coeff_const = 1.0) const
    {
      double anisotropicdiam = 0.0;
      auto vertices_index = ma->GetElVertices (ei);
      for (auto vertex1 : vertices_index)
        {
          for (auto vertex2 : vertices_index)
            {
              Vec<D> v = ma->GetPoint<D> (vertex2) - ma->GetPoint<D> (vertex1);
              vtimes (v, coeff_const);
              anisotropicdiam = max (anisotropicdiam, sqrt (L2Norm2 (v)));
            }
        }
      return anisotropicdiam * usescale + (usescale == 0);
    }
    template <int D> Vec<D> ElCenter (ElementId ei) const
    {
      Vec<D> center = 0;
      auto vertices_index = ma->GetElVertices (ei);
      for (auto vertex : vertices_index)
        center += ma->GetPoint<D> (vertex);
      center *= (1.0 / vertices_index.Size ()) * useshift;
      return center;
    }
  };

  //////////////////////////// Trefftz basis ////////////////////////////

  template <int D> class TWaveBasis : public PolBasis
  {
  public:
    TWaveBasis () { ; }
    ~TWaveBasis () { ; }
    static CSR Basis (int ord, int basistype = 0, int fowave = 0);
  };

  template <int D> class THeatBasis : public PolBasis
  {
  public:
    THeatBasis () { ; }
    ~THeatBasis () { ; }
    static CSR Basis (int ord, int basistype = 0, int fowave = 0);
  };

  template <int D> class TLapBasis : public PolBasis
  {
  public:
    TLapBasis () { ; }
    ~TLapBasis () { ; }
    static CSR Basis (int ord, int basistype = 0);
  };

  template <int D> class FOTWaveBasis : public PolBasis
  {
  public:
    FOTWaveBasis () { ; }
    ~FOTWaveBasis () { ; }
    static CSR Basis (int ord, int rdim);
  };

  //////////////////////////// quasi-Trefftz basis ////////////////////////////

  template <int D> class QTEllipticBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    std::map<string, CSR> gtbstore;

    Vector<shared_ptr<CoefficientFunction>> AAder;
    Vector<shared_ptr<CoefficientFunction>> BBder;
    Vector<shared_ptr<CoefficientFunction>> CCder;
    Vector<shared_ptr<CoefficientFunction>> FFder;

  public:
    QTEllipticBasis (int aorder, shared_ptr<CoefficientFunction> coeffA,
                     shared_ptr<CoefficientFunction> coeffB,
                     shared_ptr<CoefficientFunction> coeffC)
        : PolBasis (aorder)
    {
      if (!coeffA)
        coeffA = make_shared<ConstantCoefficientFunction> (1);
      if (!coeffB)
        coeffB = make_shared<ConstantCoefficientFunction> (0);
      if (!coeffC)
        coeffC = make_shared<ConstantCoefficientFunction> (0);

      this->ComputeDerivs<D> (order - 1, coeffA, AAder);
      this->ComputeDerivs<D> (order - 1, coeffB, BBder);
      this->ComputeDerivs<D> (order - 1, coeffC, CCder);
    }
    ~QTEllipticBasis () { ; }
    CSR Basis (Vec<D> ElCenter, double elsize = 1.0);
    void SetRHS (shared_ptr<CoefficientFunction> coeffF) override
    {
      this->ComputeDerivs<D> (order, coeffF, FFder);
    }
    void GetParticularSolution (Vec<D> ElCenter, Vec<D> elsize,
                                FlatVector<> sol, LocalHeap &lh) override;
  };

  template <int D> class QTWaveBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    std::map<string, CSR> gtbstore;
    Vector<shared_ptr<CoefficientFunction>> AAder;
    Vector<shared_ptr<CoefficientFunction>> BBder;

  public:
    QTWaveBasis () { ; }

    QTWaveBasis (int aorder, shared_ptr<CoefficientFunction> coeffA,
                 shared_ptr<CoefficientFunction> coeffB)
        : PolBasis (aorder)
    {
      if (!coeffA)
        coeffA = make_shared<ConstantCoefficientFunction> (1);
      if (!coeffB)
        coeffB = make_shared<ConstantCoefficientFunction> (1);

      // if (i == 0 && (eqtyp == "qtwave" || eqtyp == "foqtwave"))
      //  coeff = UnaryOpCF(coeffB/coeffA,GenericSqrt());
      shared_ptr<CoefficientFunction> coeffAA
          = make_shared<ConstantCoefficientFunction> (1) / (coeffA * coeffA);
      // shared_ptr<CoefficientFunction> coeffx =
      // make_shared<ConstantCoefficientFunction> (1)
      /// (coeffA * coeffA);

      this->ComputeDerivs<D - 1> (order - 2, coeffAA, AAder);
      this->ComputeDerivs<D - 1> (order - 1, coeffB, BBder);
    }
    ~QTWaveBasis () { ; }

    CSR
    Basis (int ord, Vec<D> ElCenter, double elsize = 1.0, int basistype = 0);
  };

  template <int D> class FOQTWaveBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    Vec<D, std::map<string, CSR>> gtbstore;
    Vector<shared_ptr<CoefficientFunction>> AAder;
    Vector<shared_ptr<CoefficientFunction>> BBder;

  public:
    FOQTWaveBasis () { ; }

    FOQTWaveBasis (int aorder, shared_ptr<CoefficientFunction> coeffA,
                   shared_ptr<CoefficientFunction> coeffB)
        : PolBasis (aorder)
    {
      if (!coeffA)
        coeffA = make_shared<ConstantCoefficientFunction> (1);
      if (!coeffB)
        coeffB = make_shared<ConstantCoefficientFunction> (1);

      shared_ptr<CoefficientFunction> coeffAA
          = make_shared<ConstantCoefficientFunction> (1) / (coeffA * coeffA);

      this->ComputeDerivs<D - 1> (order - 1, coeffAA, AAder);
      this->ComputeDerivs<D - 1> (order - 1, coeffB, BBder);
    }
    ~FOQTWaveBasis () { ; }

    CSR Basis (int ord, int rdim, Vec<D> ElCenter, double elsize = 1.0);
  };

  template <int D> class QTHeatBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    std::map<string, CSR> gtbstore;

    Vector<shared_ptr<CoefficientFunction>> AAder;
    Vector<shared_ptr<CoefficientFunction>> FFder;

  public:
    QTHeatBasis (int aorder, shared_ptr<CoefficientFunction> coeffA)
        : PolBasis (aorder)
    {
      if (!coeffA)
        coeffA = make_shared<ConstantCoefficientFunction> (1);

      this->ComputeDerivs<D> (order - 1, coeffA, AAder);
    }
    ~QTHeatBasis () { ; }
    CSR Basis (Vec<D> ElCenter, double hx, double ht);
    void SetRHS (shared_ptr<CoefficientFunction> coeffF) override
    {
      this->ComputeDerivs<D> (order, coeffF, FFder);
    }
    void GetParticularSolution (Vec<D> ElCenter, Vec<D> elsize,
                                FlatVector<> sol, LocalHeap &lh) override;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
