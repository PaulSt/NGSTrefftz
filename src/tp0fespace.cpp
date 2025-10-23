#include <comp.hpp>
#include <python_comp.hpp>
#include "tp0fespace.hpp"

namespace ngcomp
{

  void
  TP0FE ::CalcShape (const IntegrationPoint &ip, BareSliceVector<> shape) const
  {
    double x = ip.Point ()[0], y = ip.Point ()[1];
    double xi = (2 * x - 1);
    double eta = (2 * y - 1);

    if (zero_axis == 1)
      swap (xi, eta);

    size_t p = order;
    size_t q = order - 2;

    STACK_ARRAY (double, mem, p + q + 2);
    double *polx = &mem[0];
    double *poly = &mem[p + 1];

    LegendrePolynomial (p, xi, polx);
    LegendrePolynomial (q, eta, poly);

    for (size_t i = 0, ii = 0; i <= p; i++)
      for (size_t j = 0; j <= q; j++)
        shape[ii++] = polx[i] * poly[j] * (1 - eta) * (1 + eta);
  }

  void TP0FE ::CalcDShape (const IntegrationPoint &ip,
                           BareSliceMatrix<> dshape) const
  {
    AutoDiff<2> x (ip.Point ()[0], 0);
    AutoDiff<2> y (ip.Point ()[1], 1);
    AutoDiff<2> xi = (2 * x - 1);
    AutoDiff<2> eta = (2 * y - 1);

    if (zero_axis == 1)
      swap (xi, eta);

    size_t p = order;
    size_t q = order - 2;

    STACK_ARRAY (AutoDiff<2>, mem, p + q + 2);
    AutoDiff<2> *polx = &mem[0];
    AutoDiff<2> *poly = &mem[p + 1];

    LegendrePolynomial (p, xi, polx);
    LegendrePolynomial (q, eta, poly);

    for (size_t i = 0, ii = 0; i <= p; i++)
      for (size_t j = 0; j <= q; j++)
        {
          AutoDiff<2> shape = polx[i] * poly[j] * (1 - eta) * (1 + eta);
          dshape (ii, 0) = shape.DValue (0);
          dshape (ii++, 1) = shape.DValue (1);
        }
  }

  TP0FESpace ::TP0FESpace (shared_ptr<MeshAccess> ama, const Flags &flags)
      : FESpace (ama, flags)
  {
    type = "TP0FESpace";

    order = int (flags.GetNumFlag ("order", 3));

    if (ma->GetDimension () == 2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<MyDiffOpId>> ();
        flux_evaluator[VOL]
            = make_shared<T_DifferentialOperator<MyDiffOpGradient>> ();
      }
    else
      {
        throw Exception ("TP0FESpace implemented only in 2D");
      }
  }

  DocInfo TP0FESpace ::GetDocu ()
  {
    auto docu = FESpace::GetDocu ();
    return docu;
  }

  int TP0FESpace ::LocalNDof (ELEMENT_TYPE eltype, int order) const
  {
    switch (eltype)
      {
      case ET_TRIG:
        return (order - 1) * order / 2;
      case ET_QUAD:
        return (order + 1) * (order - 1);
      default:
        throw Exception ("TP0FESpace::LocalNDof: element type not supported");
      }
  }

  void TP0FESpace ::Update ()
  {
    FESpace::Update ();
    if (order_policy == OLDSTYLE_ORDER)
      order_policy = CONSTANT_ORDER;

    nel = ma->GetNE ();
    first_element_dof.SetSize (nel);
    if (order_policy == CONSTANT_ORDER || order_inner.Size () == 0)
      {
        order_inner.SetSize0 ();
        ndof = 0;
        for (int i = 0; i < nel; i++)
          {
            first_element_dof[i] = ndof;
            ndof += LocalNDof (ma->GetElement (ElementId (VOL, i)).GetType (),
                               this->order);
          }
        // ndof = LocalNDof (this->order) * nel;
        SetNDof (ndof);
      }
    else if (order_policy == VARIABLE_ORDER)
      {
        ndof = 0;
        for (int i = 0; i < nel; i++)
          {
            first_element_dof[i] = ndof;
            ndof += LocalNDof (ma->GetElement (ElementId (VOL, i)).GetType (),
                               order_inner[i]);
          }
        SetNDof (ndof);
      }
    else
      throw Exception ("TP0FESpace: invalid order policy");

    UpdateCouplingDofArray ();
  }

  void TP0FESpace ::SetOrder (NodeId ni, int norder)
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception ("In TP0FESpace::SetOrder. Order policy is "
                       "constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;

    if (order < 0)
      order = 0;

    if (CoDimension (ni.GetType (), ma->GetDimension ()) == 0)
      {
        if (order_inner.Size () == 0)
          {
            order_inner.SetSize (ma->GetNE ());
            order_inner = this->order;
          }
        order_inner[ni.GetNr ()] = norder;
      }
    else
      throw Exception (
          "TP0FESpace::SetOrder requires NodeType of codimension 0!");
  }

  int TP0FESpace ::GetOrder (NodeId ni) const
  {
    if (CoDimension (ni.GetType (), ma->GetDimension ()) == 0
        && ni.GetNr () < order_inner.Size ())
      return order_inner[ni.GetNr ()];
    return 0;
  }

  void TP0FESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize (0);
    if (!DefinedOn (ei) || ei.VB () != VOL)
      return;
    int first_dof = first_element_dof[ei.Nr ()];
    int local_ndof;
    if (order_inner.Size () > 0)
      local_ndof
          = LocalNDof (ma->GetElement (ei).GetType (), order_inner[ei.Nr ()]);
    else
      local_ndof = LocalNDof (ma->GetElement (ei).GetType (), this->order);
    for (int j = first_dof; j < first_dof + local_ndof; j++)
      dnums.Append (j);
  }

  void TP0FESpace ::UpdateCouplingDofArray ()
  {
    ctofdof.SetSize (ndof);
    for (auto i : Range (ma->GetNE ()))
      {
        bool definedon = DefinedOn (ElementId (VOL, i));
        Array<DofId> dofs;
        GetDofNrs (i, dofs);
        for (auto r : dofs)
          ctofdof[r] = definedon ? LOCAL_DOF : UNUSED_DOF;
      }
  }

  FiniteElement &TP0FESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement (ei);
    ELEMENT_TYPE eltype = ngel.GetType ();
    int D = ma->GetDimension ();
    int order = order_inner.Size () > 0 ? order_inner[ei.Nr ()] : this->order;

    if (ei.IsVolume () && order > 1)
      {
        if (D == 2)
          switch (eltype)
            {
            case ET_QUAD:
              {
                ElementTransformation &trafo = ma->GetTrafo (ei, alloc);
                IntegrationPoint ip;
                Mat<2, 2> jac;
                trafo.CalcJacobian (ip, jac);
                double xscale = L2Norm (jac.Col (0));
                double yscale = L2Norm (jac.Col (1));
                int zero_axis = yscale > xscale ? 1 : 0;

                return *(new (alloc)
                             TP0FE (LocalNDof (ET_QUAD, order), order,
                                    IVec<4> (ngel.Vertices ()), zero_axis));
                break;
              }
            case ET_TRIG:
              {
                Ngs_Element ngel
                    = ma->GetElement<ET_trait<ET_TRIG>::DIM, VOL> (ei.Nr ());
                L2HighOrderFE<ET_TRIG> *hofe
                    = new (alloc) L2HighOrderFE<ET_TRIG> ();
                hofe->SetVertexNumbers (ngel.vertices);
                hofe->L2HighOrderFE<ET_TRIG>::SetOrder (order - 2);
                hofe->L2HighOrderFE<ET_TRIG>::ComputeNDof ();
                return *hofe;
                break;
              }
            default:
              throw Exception ("eltype not supported in TP0FE");
            }
        else
          throw Exception ("dimension not supported in TP0FE");
      }
    try
      {
        return SwitchET<ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD> (
            eltype, [&alloc] (auto et) -> FiniteElement & {
              return *new (alloc) DummyFE<et.ElementType ()>;
            });
      }
    catch (Exception &e)
      {
        throw Exception ("illegal element type in Trefftz::GetSurfaceFE");
      }
  }

  static RegisterFESpace<TP0FESpace> initifes ("TP0FESpace");
}

#ifdef NGS_PYTHON
void ExportTP0FESpace (py::module m)
{
  using namespace ngcomp;

  ExportFESpace<TP0FESpace> (m, "TP0FESpace");
}
#endif // NGS_PYTHON
