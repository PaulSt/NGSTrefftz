#ifndef FILE_TP0FESPACE_HPP
#define FILE_TP0FESPACE_HPP

#include <fem.hpp>

namespace ngcomp
{
  class TP0FE : public FiniteElement, public ET_trait<ET_QUAD>
  {
    FlatArray<int> vnums;
    int zero_axis;

  public:
    TP0FE (int ndof, int order, FlatArray<int> avnums, int azero_axis)
        : FiniteElement (ndof, order), vnums (avnums), zero_axis (azero_axis)
    {
      ;
    }
    ELEMENT_TYPE ElementType () const override { return ET_QUAD; }

    void CalcShape (const IntegrationPoint &ip, BareSliceVector<> shape) const;
    void
    CalcDShape (const IntegrationPoint &ip, BareSliceMatrix<> dshape) const;
  };

  class MyDiffOpId : public DiffOp<MyDiffOpId>
  {
  public:
    static constexpr int DIM = 1;
    static constexpr int DIM_SPACE = 2;
    static constexpr int DIM_ELEMENT = 2;
    static constexpr int DIM_DMAT = 1;
    static constexpr int DIFFORDER = 0;

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement &fel, const MIP &mip,
                                MAT &mat, LocalHeap &lh)
    {
      HeapReset hr (lh);
      // Cast (fel).CalcShape (mip.IP (), mat.Row (0));
      switch (fel.ElementType ())
        {
        case ET_QUAD:
          static_cast<const TP0FE &> (fel).CalcShape (mip.IP (), mat.Row (0));
          break;
        case ET_TRIG:
          static_cast<const L2HighOrderFE<ET_TRIG> &> (fel).CalcShape (
              mip.IP (), mat.Row (0));
          break;
        default:
          throw Exception ("MyDiffOpId: " + ToString (fel.ElementType ())
                           + " not supported");
        }
    }
  };

  class MyDiffOpGradient : public DiffOp<MyDiffOpGradient>
  {
  public:
    static constexpr int DIM = 1;
    static constexpr int DIM_SPACE = 2;
    static constexpr int DIM_ELEMENT = 2;
    static constexpr int DIM_DMAT = 2;
    static constexpr int DIFFORDER = 1;

    static string Name () { return "grad"; }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement &fel, const MIP &mip,
                                MAT &mat, LocalHeap &lh)
    {
      HeapReset hr (lh);
      FlatMatrixFixWidth<2> dshape (fel.GetNDof (), lh);
      Cast (fel).CalcDShape (mip.IP (), dshape);
      mat = Trans (dshape * mip.GetJacobianInverse ());
    }

    static const TP0FE &Cast (const FiniteElement &fel)
    {
      return static_cast<const TP0FE &> (fel);
    }
  };

  class TP0FESpace : public FESpace
  {
    int nel;
    Array<int> order_inner;
    Array<int> first_element_dof;
    bool allow_both_axes_zero = false;

    int LocalNDof (ElementId el, int order) const;

    int GetZeroAxis (ElementId ei) const;

  public:
    TP0FESpace (shared_ptr<MeshAccess> ama, const Flags &flags);

    string GetClassName () const override { return "TP0FESpace"; }

    static DocInfo GetDocu ();

    void Update () override;
    virtual void UpdateCouplingDofArray () override;

    virtual void SetOrder (NodeId ni, int norder) override;
    virtual int GetOrder (NodeId ni) const override;

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTP0FESpace (py::module m);
#endif // NGS_PYTHON

#endif
