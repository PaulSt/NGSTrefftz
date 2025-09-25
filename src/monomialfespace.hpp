#ifndef FILE_MONOMIALFESPACE_HPP
#define FILE_MONOMIALFESPACE_HPP

#include "scalarmappedfe.hpp"

#include <comp.hpp>

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

    void Update () override;
    virtual void UpdateCouplingDofArray () override;

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;

    string GetClassName () const override { return "monomialfespace"; }
    static DocInfo GetDocu ();

  protected:
    template <int D> CSR MonomialBasis (int ord) const
    {
      CSR tb;
      const int npoly = BinCoeff (D + ord, ord);
      Matrix<> basis (npoly, npoly);
      basis = 0.0;
      for (int i = 0; i < npoly; i++)
        basis (i, i) = 1.0;

      MatToCSR (basis, tb);
      return tb;
    }
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
      return anisotropicdiam;
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
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMonomialFESpace (py::module m);
#endif // NGS_PYTHON

#endif
