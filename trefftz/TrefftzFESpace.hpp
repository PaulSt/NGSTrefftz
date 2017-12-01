#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP

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

  public:
    /*
      constructor.
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags);

    virtual string GetClassName () const { return "TrefftzFESpace"; }

    virtual void Update (LocalHeap &lh);

    virtual void GetDofNrs (ElementId ei, Array<DofId> &dnums) const;
    virtual FiniteElement &GetFE (ElementId ei, Allocator &alloc) const;

    // int GetDim() const { return local_ndof; }

    template <int D> double Adiam (ElementId ei) const;

    virtual size_t GetNDof () const throw () override { return ndof; }

  protected:
    template <ELEMENT_TYPE ET>
    FiniteElement &T_GetFE (int elnr, Allocator &alloc) const;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
