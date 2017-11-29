#ifndef FILE_L2HOFESPACEE
#define FILE_L2HOFESPACEE

/*********************************************************************/
/* File:   l2hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   23.Feb. 2003                                              */
/*********************************************************************/




namespace ngcomp
{

  /**
     High Order Finite Element Space for L2 (element by element)
  */

  class NGS_DLL_HEADER L2HighOrderFESpace2 : public FESpace
  {
  protected:
    // Number of Elements
    int nel;
    // Degrees of Freedom
    int ndof;
    // Levels
    Array<int> ndlevel;
    // if order is relative to mesh order
    bool var_order;
    // variable order is set to mesh_order + rel_order
    int rel_order;
    // order of elements
    Array<INT<3> > order_inner;
    // table of first element dofnumber
    Array<DofId> first_element_dof;
    bool all_dofs_together;
  public:

    L2HighOrderFESpace2 (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~L2HighOrderFESpace2 ();
    // Create if order=0 ElementFESpace Constructor, else L2HOFE
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);
    // Creates also for order=0 a L2HighOrderFESpace2
    static shared_ptr<FESpace> CreateHO (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<L2HighOrderFESpace2> (ma, flags);
    }

    virtual string GetClassName () const override
    {
      return "L2HighOrderFESpace2";
    }

    bool AllDofsTogether(){return all_dofs_together;};
    ///
    virtual void Update(LocalHeap & lh) override;
    ///
    virtual void UpdateDofTables();
    ///
    // virtual void UpdateCouplingDofArray();
    ///
    virtual size_t GetNDof () const throw() override;
    ///
    virtual size_t GetNDofLevel (int level) const override;
    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    using FESpace::GetFE;
    // virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const override;
    // ///
    // virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const override;
    ///
    // virtual const FiniteElement & GetFacetFE (int fnr, LocalHeap & lh) const;

    // virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///
    // virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
    ///

    // virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    // virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    // virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    // virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;


    auto GetElementDofs (size_t nr) const
    {
      return Range (first_element_dof[nr], first_element_dof[nr+1]);
    }

    // virtual void SolveM (CoefficientFunction & rho, BaseVector & vec,
    //                      LocalHeap & lh) const override;


  protected:

    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE (int elnr, Allocator & alloc) const;
  };

}

#endif
