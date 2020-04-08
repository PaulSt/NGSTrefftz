#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP
#include "trefftzwavefe.hpp"

namespace ngcomp
{

  class TrefftzFESpace : public FESpace
  {
    int D;
    int fullD;
    int order;
    size_t ndof;
    int nel;
    int nvert;
    int local_ndof;
    float c = 1;
    int useshift = 1;
    int basistype;
    Array<double> gamma;
    int gppword;

  public:
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags);

    string GetClassName () const override { return "trefftz"; }

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;

    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;

    size_t GetNDof () const override { return ndof; }

    static DocInfo GetDocu ();

  protected:
    template <int D> double Adiam (ElementId ei) const;

    template <int D> Vec<D + 1> ElCenter (ElementId ei) const;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
