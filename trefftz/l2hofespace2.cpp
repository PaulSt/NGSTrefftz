/*********************************************************************/
/* File:   l2hofespace.cpp                                         */
/* Author: Start                                                     */
/* Date:   24. Feb. 2003                                             */
/*********************************************************************/

/**
   High Order Finite Element Space for L2
*/

/* ***********************************************
To do: *Internal External Dofs (eliminate internal)
       *Flag for low-order dofs eliminated ...
************************* */

#include <comp.hpp>
#include <multigrid.hpp>
#include "l2hofespace2.hpp"

using namespace ngmg;

namespace ngcomp
{

  // class BlockDifferentialOperatorId : public BlockDifferentialOperator
  // {
  // public:
  //   using BlockDifferentialOperator::BlockDifferentialOperator;
  //
  //   virtual void Apply (const FiniteElement & fel,
  //                       const SIMD_BaseMappedIntegrationRule & mir,
  //                       BareSliceVector<double> x,
  //                       BareSliceMatrix<SIMD<double>> flux) const
  //   {
  //     if (comp == -1)
  //       static_cast<const BaseScalarFiniteElement&> (fel).
  //         Evaluate(mir.IR(), SliceMatrix<double> (fel.GetNDof(), dim, dim,
  //         &x(0)), flux);
  //     else
  //       diffop->Apply(fel, mir, x.Slice(comp, dim),
  //       flux.RowSlice(comp,dim));
  //   }
  //   virtual void
  //   AddTrans (const FiniteElement & fel,
  //             const SIMD_BaseMappedIntegrationRule & mir,
  //             BareSliceMatrix<SIMD<double>> flux,
  //             BareSliceVector<double> x) const
  //   {
  //   if (comp == -1)
  //     static_cast<const BaseScalarFiniteElement&> (fel).
  //       AddTrans(mir.IR(), flux, SliceMatrix<double> (fel.GetNDof(), dim,
  //       dim, &x(0)));
  //   else
  //     diffop->AddTrans(fel, mir, flux.RowSlice(comp,dim),
  //     x.Slice(comp,dim));
  //   }
  //
  // };

  L2HighOrderFESpace2 ::L2HighOrderFESpace2 (shared_ptr<MeshAccess> ama,
                                             const Flags &flags,
                                             bool parseflags)
      : FESpace (ama, flags)
  {
    name = "L2HighOrderFESpace2(l2ho)";

    // defined flags
    DefineNumFlag ("relorder");
    DefineDefineFlag ("l2ho");
    DefineDefineFlag ("all_dofs_together");

    if (parseflags)
      CheckFlags (flags);

    var_order = 0;
    order = int (flags.GetNumFlag ("order", 0));

    switch (ma->GetDimension ())
      {
      case 1:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<1>>> ();
          flux_evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpGradient<1>>> ();
          // evaluator[BND] =
          // make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>();
          break;
        }
      case 2:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>> ();
          flux_evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpGradient<2>>> ();
          // evaluator[BND] =
          // make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
          break;
        }
      case 3:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>> ();
          flux_evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpGradient<3>>> ();
          // evaluator[BND] =
          // make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
          break;
        }
      }

    all_dofs_together = true;
  }

  L2HighOrderFESpace2 ::~L2HighOrderFESpace2 () { ; }

  shared_ptr<FESpace>
  L2HighOrderFESpace2 ::Create (shared_ptr<MeshAccess> ma, const Flags &flags)
  {
    int order = int (flags.GetNumFlag ("order", 0));
    if (order == 0)
      return make_shared<ElementFESpace> (ma, flags);
    else
      return make_shared<L2HighOrderFESpace2> (ma, flags, true);
  }

  void L2HighOrderFESpace2 ::Update (LocalHeap &lh)
  {
    FESpace::Update (lh);
    if (low_order_space)
      low_order_space->Update (lh);

    nel = ma->GetNE ();
    order_inner.SetSize (nel);

    order_inner = INT<3> (order);

    if (var_order)
      for (int i = 0; i < nel; i++)
        order_inner[i] = ma->GetElOrders (i) + INT<3> (rel_order);

    for (int i = 0; i < nel; i++)
      {
        ElementId ei (VOL, i);
        order_inner[i]
            = order_inner[i] + INT<3> (et_bonus_order[ma->GetElType (ei)]);
        order_inner[i] = Max (order_inner[i], INT<3> (0));
        if (!DefinedOn (VOL, ma->GetElIndex (ei)))
          order_inner[i] = 0;
      }
    if (print)
      *testout << " order_inner (l2ho) " << order_inner << endl;

    UpdateDofTables ();
  }

  void L2HighOrderFESpace2 ::UpdateDofTables ()
  {
    ndof = all_dofs_together ? 0 : nel;
    first_element_dof.SetSize (nel + 1);
    for (int i = 0; i < nel; i++)
      {
        ElementId ei (VOL, i);
        first_element_dof[i] = ndof;
        INT<3> pi = order_inner[i];
        switch (ma->GetElType (ei))
          {
          case ET_SEGM:
            ndof += pi[0] + 1;
            break;
          case ET_TRIG:
            ndof += (pi[0] + 1) * (pi[0] + 2) / 2;
            break;
          case ET_QUAD:
            ndof += (pi[0] + 1) * (pi[1] + 1);
            break;
          case ET_TET:
            ndof += (pi[0] + 1) * (pi[0] + 2) * (pi[0] + 3) / 6;
            break;
          case ET_PRISM:
            ndof += (pi[0] + 1) * (pi[0] + 2) * (pi[2] + 1) / 2;
            break;
          case ET_PYRAMID:
            ndof += 5 + 8 * (pi[0] - 1) + 2 * (pi[0] - 1) * (pi[0] - 2)
                    + (pi[0] - 1) * (pi[0] - 1)
                    + (pi[0] - 1) * (pi[0] - 2) * (2 * pi[0] - 3) / 6;
            break;
          case ET_HEX:
            ndof += (pi[0] + 1) * (pi[1] + 1) * (pi[2] + 1);
            break;
          default: // for the compiler
            break;
          }
        if (!all_dofs_together)
          ndof--; // subtract constant
      }
    first_element_dof[nel] = ndof;
    // cout << first_element_dof << endl;

    if (print)
      *testout << " first_element dof (l2hofe) " << first_element_dof << endl;
    cout << "===========================update ndof: " << ndof
         << " ndof/el: " << ndof / nel << endl;
  }

  FiniteElement &
  L2HighOrderFESpace2 ::GetFE (ElementId ei, Allocator &alloc) const
  {
    if (ei.VB () == BBND)
      throw Exception ("BBND not available in L2HighOrderFESpace2");
    if (ei.IsVolume ())
      {
        int elnr = ei.Nr ();
        Ngs_Element ngel = ma->GetElement (ei);
        ELEMENT_TYPE eltype = ngel.GetType ();

        if (!DefinedOn (ngel))
          {
            switch (eltype)
              {
              case ET_POINT:
                return *new (alloc) ScalarDummyFE<ET_POINT> ();
                break;
              case ET_SEGM:
                return *new (alloc) ScalarDummyFE<ET_SEGM> ();
                break;
              case ET_TRIG:
                return *new (alloc) ScalarDummyFE<ET_TRIG> ();
                break;
              case ET_QUAD:
                return *new (alloc) ScalarDummyFE<ET_QUAD> ();
                break;
              case ET_TET:
                return *new (alloc) ScalarDummyFE<ET_TET> ();
                break;
              case ET_PYRAMID:
                return *new (alloc) ScalarDummyFE<ET_PYRAMID> ();
                break;
              case ET_PRISM:
                return *new (alloc) ScalarDummyFE<ET_PRISM> ();
                break;
              case ET_HEX:
                return *new (alloc) ScalarDummyFE<ET_HEX> ();
                break;
              }
            cout << "returning scalar dummy" << endl;
          }

        if (eltype == ET_TRIG)
          {
            // return *CreateL2HighOrderFE<ET_TRIG> (order, vnums, alloc);
            // cout << "volume trig" << endl;
            return *CreateL2HighOrderFE<ET_TRIG> (
                order, INT<3> (ngel.Vertices ()), alloc);
          }

        if (eltype == ET_TET)
          {
            // cout << "volume tet" << endl;
            return *CreateL2HighOrderFE<ET_TET> (
                order, INT<4> (ngel.Vertices ()), alloc);
          }

        switch (eltype)
          {
            cout << "volume switch" << endl;
          case ET_SEGM:
            return T_GetFE<ET_SEGM> (elnr, alloc);
          // case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, alloc);
          case ET_QUAD:
            return T_GetFE<ET_QUAD> (elnr, alloc);
          // case ET_TET:     return T_GetFE<ET_TET> (elnr, alloc);
          case ET_PRISM:
            return T_GetFE<ET_PRISM> (elnr, alloc);
          case ET_PYRAMID:
            return T_GetFE<ET_PYRAMID> (elnr, alloc);
          case ET_HEX:
            return T_GetFE<ET_HEX> (elnr, alloc);
          default:
            throw Exception ("illegal element in L2HoFeSpace::GetFE");
          }
      }
    else
      {
        int elnr = ei.Nr ();
        switch (ma->GetElType (ei))
          {
            cout << "else dummy" << endl;
          case ET_POINT:
            return *new (alloc) DummyFE<ET_POINT>;
          case ET_SEGM:
            return *new (alloc) DummyFE<ET_SEGM>;
            break;
          case ET_TRIG:
            return *new (alloc) DummyFE<ET_TRIG>;
            break;
          case ET_QUAD:
            return *new (alloc) DummyFE<ET_QUAD>;
            break;

          default:
            stringstream str;
            str << "FESpace " << GetClassName ()
                << ", undefined surface eltype " << ma->GetElType (ei)
                << ", order = " << order << endl;
            throw Exception (str.str ());
          }
      }
  }

  template <ELEMENT_TYPE ET>
  FiniteElement &L2HighOrderFESpace2 ::T_GetFE (int elnr, Allocator &lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM, VOL> (elnr);
    L2HighOrderFE<ET> *hofe = new (lh) L2HighOrderFE<ET> ();

    hofe->SetVertexNumbers (ngel.vertices);
    hofe->L2HighOrderFE<ET>::SetOrder (order_inner[elnr]);
    hofe->L2HighOrderFE<ET>::ComputeNDof ();

    return *hofe;
  }

  size_t L2HighOrderFESpace2 ::GetNDof () const throw ()
  {
    // cout << "ndof: " << ndof << endl;
    return ndof;
  }

  size_t L2HighOrderFESpace2 ::GetNDofLevel (int level) const
  {
    cout << "ndof level: " << ndlevel[level] << endl;
    return ndlevel[level];
  }

  void L2HighOrderFESpace2 ::GetDofNrs (ElementId ei, Array<int> &dnums) const
  {
    dnums.SetSize0 ();
    if (!DefinedOn (ei) || ei.VB () != VOL)
      return;

    auto eldofs = GetElementDofs (ei.Nr ());
    size_t size = eldofs.Size ();
    size_t base = all_dofs_together ? 0 : 1;
    size += base;
    dnums.SetSize (size);
    if (!all_dofs_together)
      dnums[0] = ei.Nr ();
    dnums.Range (base, size) = eldofs;
    // cout << "GetDofNrs: ei.Nr() = " << ei.Nr() << " ndof: " << ndof <<"
    // dnums: \n" << dnums << endl <<
    // "================================================" << endl ;
  }

  // register FESpaces
  namespace l2hofespace_cpp
  {
    class Init
    {
    public:
      Init ();
    };

    Init::Init ()
    {
      GetFESpaceClasses ().AddFESpace ("l22", L2HighOrderFESpace2::Create);
      // GetFESpaceClasses().AddFESpace ("l22ho",
      // L2HighOrderFESpace2::CreateHO); GetFESpaceClasses().AddFESpace
      // ("l22surf", L2SurfaceHighOrderFESpace::Create);
    }

    Init init;
  }
}
