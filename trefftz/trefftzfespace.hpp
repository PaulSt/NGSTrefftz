#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP
#include "trefftzelement.hpp"

namespace ngcomp
{

  class TrefftzFESpace : public FESpace
  {
    int D;
    int order;
    size_t ndof;
    int nel;
    int nvert;
    int local_ndof;
    float c = 1;
    bool testshift = false;
    int basistype;

  public:
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags);

    virtual string GetClassName () const { return "TrefftzFESpace"; }

    virtual void Update (LocalHeap &lh);

    virtual void GetDofNrs (ElementId ei, Array<DofId> &dnums) const;

    virtual FiniteElement &GetFE (ElementId ei, Allocator &alloc) const;

    virtual size_t GetNDof () const throw () override { return ndof; }

  protected:
    template <int D> double Adiam (ElementId ei) const;

    template <int D> Vec<D> ElCenter (ElementId ei) const;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
