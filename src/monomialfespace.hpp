#ifndef FILE_MONOMIALFESPACE_HPP
#define FILE_MONOMIALFESPACE_HPP

#include "scalarmappedfe.hpp"

namespace ngcomp
{

  class MonomialFESpace : public FESpace
  {
    int D;
    int order;
    int nel;
    int local_ndof;
    int useshift = 1;
    int usescale = 1;
    shared_ptr<CoefficientFunction> coeff_cf = nullptr;
    CSR basismat;

  public:
    MonomialFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
                     bool checkflags = false);

    void SetCoeff (shared_ptr<CoefficientFunction> acoeff_cf)
    {
      coeff_cf = acoeff_cf;
    }

    string GetClassName () const override { return "monomialfespace"; }

    void Update () override;

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;

    virtual void UpdateCouplingDofArray () override;

    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;

    static DocInfo GetDocu ();

  protected:
    template <int D> CSR MonomialBasis (int ord) const
    {
      CSR tb;
      const int npoly = BinCoeff (D + 1 + ord, ord);
      Matrix<> basis (npoly, npoly);
      basis = 0.0;
      for (int i = 0; i < npoly; i++)
        basis (i, i) = 1.0;

      MatToCSR (basis, tb);
      return tb;
    }

    template <int D> double Adiam (ElementId ei) const
    {
      if (usescale == 0)
        return 1.0;
      LocalHeap lh (1000 * 1000);
      double anisotropicdiam = 0.0;
      auto vertices_index = ma->GetElVertices (ei);

      for (auto vertex1 : vertices_index)
        {
          for (auto vertex2 : vertices_index)
            {
              Vec<D> v1 = ma->GetPoint<D> (vertex1);
              Vec<D> v2 = ma->GetPoint<D> (vertex2);
              IntegrationRule ir (ma->GetElType (ei), 0);
              ElementTransformation &trafo = ma->GetTrafo (ei, lh);
              MappedIntegrationPoint<D, D> mip (ir[0], trafo);
              mip.Point () = v1;
              double c1 = coeff_cf ? coeff_cf->Evaluate (mip) : 1.0;
              mip.Point () = v2;
              double c2 = coeff_cf ? coeff_cf->Evaluate (mip) : 1.0;

              anisotropicdiam = max (
                  anisotropicdiam,
                  sqrt (L2Norm2 (v1.Range (0, D - 1) - v2.Range (0, D - 1))
                        + pow (c1 * v1 (D - 1) - c2 * v2 (D - 1), 2)));
            }
        }
      return anisotropicdiam;
    }

    template <int D> Vec<D> ElCenter (ElementId ei) const
    {
      if (useshift == 0)
        return 1;
      Vec<D> center = 0;
      auto vertices_index = ma->GetElVertices (ei);
      for (auto vertex : vertices_index)
        center += ma->GetPoint<D> (vertex);
      center *= (1.0 / vertices_index.Size ());
      return center;
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMonomialFESpace (py::module m);
#endif // NGS_PYTHON

#endif
