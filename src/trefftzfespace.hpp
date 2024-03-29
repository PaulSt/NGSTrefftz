#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP

#include "scalarmappedfe.hpp"
#include "planewavefe.hpp"

namespace ngcomp
{

  class PolBasis
  {
  protected:
    int order;

  public:
    PolBasis () { ; }
    PolBasis (int aorder) : order (aorder) { ; }
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
          for (int i = 0, ii = 0; i <= order; i++)
            {
              ders (i) = coeffx;
              coeffx = coeffx->Diff (
                  MakeCoordinateCoefficientFunction (0).get (),
                  make_shared<ConstantCoefficientFunction> (1));
            }
          break;
        case 2:
          for (int i = 0, ii = 0; i <= order; i++)
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
          for (int i = 0, ii = 0; i <= order; i++)
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
    string eqtyp = "wave";
    int useshift = 1;
    int usescale = 1;
    int basistype = 0;
    shared_ptr<CoefficientFunction> coeffA = nullptr;
    shared_ptr<CoefficientFunction> coeffB = nullptr;
    shared_ptr<CoefficientFunction> coeffC = nullptr;

    CSR basismat;
    Vector<CSR> basismats;
    PolBasis *basis;

  public:
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags);
    void SetCoeff (double acoeff_const);
    void SetCoeff (shared_ptr<CoefficientFunction> acoeffA,
                   shared_ptr<CoefficientFunction> acoeffB = nullptr,
                   shared_ptr<CoefficientFunction> acoeffC = nullptr);

    shared_ptr<GridFunction>
    GetEWSolution (shared_ptr<CoefficientFunction> acoeffF);
    string GetClassName () const override { return "trefftz"; }
    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;
    size_t GetNDof () const override { return ndof; }
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
    static CSR Basis (int ord, int basistype = 0, int fowave = 0);
  };

  template <int D> class THeatBasis : public PolBasis
  {
  public:
    THeatBasis () { ; }
    static CSR Basis (int ord, int basistype = 0, int fowave = 0);
  };

  template <int D> class TLapBasis : public PolBasis
  {
  public:
    TLapBasis () { ; }
    static CSR Basis (int ord, int basistype = 0);
  };

  template <int D> class FOTWaveBasis : public PolBasis
  {
  public:
    FOTWaveBasis () { ; }
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
    CSR Basis (Vec<D> ElCenter, double elsize = 1.0);
    void SetRHS (shared_ptr<CoefficientFunction> coeffF)
    {
      this->ComputeDerivs<D> (order, coeffF, FFder);
    }
    void GetParticularSolution (Vec<D> ElCenter, double elsize,
                                FlatVector<> sol, LocalHeap &lh);
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

    CSR Basis (int ord, int rdim, Vec<D> ElCenter, double elsize = 1.0);
  };

  template <int D> class QTHeatBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    std::map<string, CSR> gtbstore;

    Vector<shared_ptr<CoefficientFunction>> AAder;

  public:
    QTHeatBasis (int aorder, shared_ptr<CoefficientFunction> coeffA)
        : PolBasis (aorder)
    {
      if (!coeffA)
        coeffA = make_shared<ConstantCoefficientFunction> (1);

      this->ComputeDerivs<D> (order - 1, coeffA, AAder);
    }
    CSR Basis (Vec<D> ElCenter, double hx, double ht);
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
