#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP


namespace ngcomp
{

  class TrefftzFESpace : public FESpace
  {
		int D;
    int order;
    int ndof;
		int nvert;
		Array<int> first_edge_dof;
    Array<int> first_cell_dof;

  public:
    /*
      constructor.
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "TrefftzFESpace";
    }

    virtual void Update(LocalHeap & lh);
    virtual size_t GetNDof () const { return ndof; }

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;

    // some new functionality our space should have in Python
    int GetNVert() { return nvert; }
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace(py::module m);
#endif // NGS_PYTHON

#endif
