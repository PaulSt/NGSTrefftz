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
    void ComputeDerivs (shared_ptr<CoefficientFunction> acoeff,
                        Matrix<shared_ptr<CoefficientFunction>> &ders);
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
    float coeff_const = 1;
    string eqtyp = "wave";
    int useshift = 1;
    int usescale = 1;
    int basistype = 0;
    shared_ptr<CoefficientFunction> coeffA = nullptr;
    shared_ptr<CoefficientFunction> coeffB = nullptr;
    shared_ptr<CoefficientFunction> coeffC = nullptr;
    Matrix<shared_ptr<CoefficientFunction>> AAder;
    Matrix<shared_ptr<CoefficientFunction>> BBder;
    Matrix<shared_ptr<CoefficientFunction>> CCder;
    CSR basismat;
    Vector<CSR> basismats;
    PolBasis *basis;

  public:
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags);
    void SetCoeff (double acoeff_const);
    void SetCoeff (shared_ptr<CoefficientFunction> acoeffA,
                   shared_ptr<CoefficientFunction> acoeffB = nullptr,
                   shared_ptr<CoefficientFunction> acoeffC = nullptr);
    string GetClassName () const override { return "trefftz"; }
    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;
    size_t GetNDof () const override { return ndof; }
    static DocInfo GetDocu ();

  protected:
    void UpdateBasis ();
    template <int D>
    double ElSize (ElementId ei, double coeff_const = 1.0) const
    {
      double anisotropicdiam = 0.0;
      auto vertices_index = ma->GetElVertices (ei);
      for (auto vertex1 : vertices_index)
        {
          for (auto vertex2 : vertices_index)
            {
              Vec<D> v1 = ma->GetPoint<D> (vertex1);
              Vec<D> v2 = ma->GetPoint<D> (vertex2);
              anisotropicdiam = max (
                  anisotropicdiam,
                  sqrt (L2Norm2 (v1.Range (0, D - 1) - v2.Range (0, D - 1))
                        + pow (coeff_const * (v1 (D - 1) - v2 (D - 1)), 2)));
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

    Matrix<shared_ptr<CoefficientFunction>> AAder;
    Matrix<shared_ptr<CoefficientFunction>> BBder;
    Matrix<shared_ptr<CoefficientFunction>> CCder;

  public:
    QTEllipticBasis (int aorder, shared_ptr<CoefficientFunction> coeffA,
                     shared_ptr<CoefficientFunction> coeffB,
                     shared_ptr<CoefficientFunction> coeffC)
        : PolBasis (aorder)
    {
      AAder.SetSize (this->order, D == 2 ? this->order : 1);
      BBder.SetSize (this->order, D == 2 ? this->order : 1);
      CCder.SetSize (this->order, D == 2 ? this->order : 1);
      if (!coeffA)
        coeffA = make_shared<ConstantCoefficientFunction> (1);
      if (!coeffB)
        coeffB = make_shared<ConstantCoefficientFunction> (0);
      if (!coeffC)
        coeffC = make_shared<ConstantCoefficientFunction> (0);

      this->ComputeDerivs (coeffA, AAder);
      this->ComputeDerivs (coeffB, BBder);
      this->ComputeDerivs (coeffC, CCder);
    }
    CSR Basis (Vec<D> ElCenter, double elsize = 1.0);
  };

  template <int D> class QTWaveBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    std::map<string, CSR> gtbstore;
    Matrix<shared_ptr<CoefficientFunction>> AAder;
    Matrix<shared_ptr<CoefficientFunction>> BBder;

  public:
    QTWaveBasis () { ; }

    QTWaveBasis (int aorder, shared_ptr<CoefficientFunction> coeffA,
                 shared_ptr<CoefficientFunction> coeffB)
        : PolBasis (aorder)
    {
      AAder.SetSize (this->order - 1, D == 2 ? this->order - 1 : 1);
      BBder.SetSize (this->order, D == 2 ? this->order : 1);
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

      this->ComputeDerivs (coeffAA, AAder);
      this->ComputeDerivs (coeffB, BBder);
    }

    CSR Basis (int ord, Vec<D + 1> ElCenter, double elsize = 1.0,
               int basistype = 0);
  };

  template <int D> class FOQTWaveBasis : public PolBasis
  {
    mutex gentrefftzbasis;
    Vec<D + 1, std::map<string, CSR>> gtbstore;
    Matrix<shared_ptr<CoefficientFunction>> AAder;
    Matrix<shared_ptr<CoefficientFunction>> BBder;

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
      AAder.SetSize (this->order - 1, D == 2 ? this->order - 1 : 1);
      BBder.SetSize (this->order, D == 2 ? this->order : 1);

      shared_ptr<CoefficientFunction> coeffAA
          = make_shared<ConstantCoefficientFunction> (1) / (coeffA * coeffA);

      this->ComputeDerivs (coeffAA, AAder);
      this->ComputeDerivs (coeffB, BBder);
    }

    CSR Basis (int ord, int rdim, Vec<D + 1> ElCenter, double elsize = 1.0);
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
