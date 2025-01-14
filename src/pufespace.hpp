#ifndef FILE_PUFESPACE_HPP
#define FILE_PUFESPACE_HPP

#include "scalarmappedfe.hpp"

#include <comp.hpp>

namespace ngcomp
{

  class PUFESpace : public FESpace
  {
    int D;
    int order;
    size_t local_ndof;
    int useshift = 1;
    int usescale = 1;
    shared_ptr<CoefficientFunction> coeff_cf = nullptr;
    CSR basismat;

  public:
    PUFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
               bool checkflags = false);

    string GetClassName () const override { return "pufespace"; }

    void Update () override;

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;

    virtual void UpdateCouplingDofArray () override;

    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;

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

    template <int D> Vec<D + 1> Adiam (ElementId ei) const
    {
      if (usescale == 0)
        return 1.0;

      Vec<D + 1> diams = 0.0;

      auto vertices = ma->GetElVertices (ei);
      for (int v = 0; v < D + 1; v++)
        {
          auto elements = ma->GetVertexElements (vertices[v]);
          for (auto element1 : elements)
            for (auto element2 : elements)
              {
                auto vertices_index1 = ma->GetElVertices (element1);
                auto vertices_index2 = ma->GetElVertices (element2);
                for (auto vertex1 : vertices_index1)
                  for (auto vertex2 : vertices_index2)
                    {
                      Vec<D> v1 = ma->GetPoint<D> (vertex1);
                      Vec<D> v2 = ma->GetPoint<D> (vertex2);
                      diams[v] = max (diams[v], L2Norm (v1 - v2));
                    }
              }
        }
      return diams;
    }

    template <int D> Vec<D + 1, Vec<D>> ElVertices (ElementId ei) const
    {
      Vec<D + 1, Vec<D>> vertices;
      auto vertices_index = ma->GetElVertices (ei);
      int i = 0;
      for (auto vertex : vertices_index)
        vertices[i++] = ma->GetPoint<D> (vertex);
      // vertices[i++] = 0.0;
      return vertices;
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportPUFESpace (py::module m);
#endif // NGS_PYTHON

#endif
