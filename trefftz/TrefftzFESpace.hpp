#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP

namespace ngcomp
{

  class TrefftzFESpace : public FESpace
  {
    int D;
    int order;
    size_t ndof;
    int nvert;
    int local_ndof;
    Array<int> first_edge_dof;
    Array<int> first_cell_dof;
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
    virtual size_t GetNDof () const { return ndof; }

    virtual void GetDofNrs (ElementId ei, Array<DofId> &dnums) const;
    virtual FiniteElement &GetFE (ElementId ei, Allocator &alloc) const;

    int GetDim () const { return local_ndof; }
    // size_t GetNDof() const {return ndof;}

    template <int D> double Adiam (ElementId ei) const;

    virtual const FiniteElement &GetFacetFE (int fnr, LocalHeap &lh) const
    {
      cout << "hello!!!!!!!!!";
    }
    virtual void GetDofRanges (ElementId ei, Array<IntRange> &dranges) const
    {
      cout << "hello!!!!!!!!!";
    }
    virtual shared_ptr<Table<int>>
    CreateSmoothingBlocks (const Flags &precflags) const override
    {
      cout << "hello!!!!!!!!!";
    }
    virtual void GetVertexDofNrs (int vnr, Array<DofId> &dnums) const override
    {
      cout << "hello!!!!!!!!!";
    }
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> &dnums) const override
    {
      cout << "hello!!!!!!!!!";
    }
    virtual void GetFaceDofNrs (int fanr, Array<DofId> &dnums) const override
    {
      cout << "hello!!!!!!!!!";
    }
    virtual void GetInnerDofNrs (int elnr, Array<DofId> &dnums) const override
    {
      cout << "hello!!!!!!!!!";
    }
    auto GetElementDofs (size_t nr) const { cout << "hello!!!!!!!!!"; }
    virtual void SolveM (CoefficientFunction &rho, BaseVector &vec,
                         LocalHeap &lh) const override
    {
      cout << "hello!!!!!!!!!";
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
