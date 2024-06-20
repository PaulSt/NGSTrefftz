#include <comp.hpp> // provides FESpace, ...
#include <python_comp.hpp>

#include "trefftzfespace.hpp"
#include "monomialfespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{
  TrefftzFESpace ::TrefftzFESpace (shared_ptr<MeshAccess> ama,
                                   const Flags &flags)
      : FESpace (ama, flags)
  {
    type = "trefftzfespace";

    D = ma->GetDimension ();

    this->dgjumps = true;
    // coeff_const = flags.GetNumFlag ("wavespeed", 1);
    // coeffA = nullptr;
    basistype = flags.GetNumFlag ("basistype", 0);
    useshift = flags.GetNumFlag ("useshift", 1);
    usescale = flags.GetNumFlag ("usescale", 1);
    DefineNumListFlag ("eq");
    eqtyp = flags.GetStringFlag ("eq");

    if (eqtyp == "fowave" || eqtyp == "foqtwave")
      local_ndof = (D)*BinCoeff (D - 1 + order, D - 1);
    else if (eqtyp == "wave" || eqtyp == "fowave_reduced")
      local_ndof = BinCoeff (D - 1 + order, order)
                   + BinCoeff (D - 1 + order - 1, order - 1)
                   - (eqtyp == "fowave_reduced");
    else if (eqtyp == "heat")
      local_ndof = BinCoeff (D - 1 + order, order);
    else if (eqtyp == "qtheat")
      local_ndof = BinCoeff (D - 1 + order, order)
                   + BinCoeff (D - 1 + order - 1, order - 1);
    else if (eqtyp == "laplace")
      local_ndof = BinCoeff (D - 1 + order, order)
                   + BinCoeff (D - 1 + order - 1, order - 1);
    else if (eqtyp == "qtelliptic")
      local_ndof = BinCoeff (D - 1 + order, order)
                   + BinCoeff (D - 1 + order - 1, order - 1);
    else if (eqtyp == "helmholtz" || eqtyp == "helmholtzconj")
      local_ndof = 2 * order + 1;
    else
      local_ndof = BinCoeff (D - 1 + order, order)
                   + BinCoeff (D - 1 + order - 1, order - 1);

    nel = ma->GetNE ();
    ndof = local_ndof * nel;

    SetDefinedOn (BND, BitArray (ma->GetNRegions (BND)).Clear ());

    // evaluators
    switch (D)
      {
      case 2:
        {
          if (eqtyp == "fowave" || eqtyp == "foqtwave")
            {
              evaluator[VOL] = make_shared<T_DifferentialOperator<
                  DiffOpMappedGradient<2, BlockMappedElement<2>>>> ();
              // flux_evaluator[VOL] =
              // make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>>();
              additional_evaluators.Set (
                  "grad",
                  make_shared<
                      T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
            }
          else if (eqtyp == "helmholtz" || eqtyp == "helmholtzconj")
            {
              evaluator[VOL] = make_shared<
                  T_DifferentialOperatorC<DiffOpMappedComplex<2>>> ();
              flux_evaluator[VOL] = make_shared<
                  T_DifferentialOperatorC<DiffOpMappedGradientComplex<2>>> ();
            }
          else
            {
              evaluator[VOL]
                  = make_shared<T_DifferentialOperator<DiffOpMapped<2>>> ();
              flux_evaluator[VOL] = make_shared<
                  T_DifferentialOperator<DiffOpMappedGradient<2>>> ();
              additional_evaluators.Set (
                  "hesse",
                  make_shared<
                      T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
            }
          break;
        }
      case 3:
        {
          if (eqtyp == "fowave" || eqtyp == "foqtwave")
            {
              evaluator[VOL] = make_shared<T_DifferentialOperator<
                  DiffOpMappedGradient<3, BlockMappedElement<3>>>> ();
              // flux_evaluator[VOL] =
              // make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>>();
              additional_evaluators.Set (
                  "grad",
                  make_shared<
                      T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
            }
          else
            {
              evaluator[VOL]
                  = make_shared<T_DifferentialOperator<DiffOpMapped<3>>> ();
              flux_evaluator[VOL] = make_shared<
                  T_DifferentialOperator<DiffOpMappedGradient<3>>> ();
              additional_evaluators.Set (
                  "hesse",
                  make_shared<
                      T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
            }
          break;
        }
      }

    UpdateBasis ();
  }

  void TrefftzFESpace ::UpdateBasis ()
  {
    string eqtyp = this->eqtyp;
    switch (D)
      {
      case 2:
        {
          if (eqtyp == "laplace")
            basismat = TLapBasis<2>::Basis (order, basistype);
          else if (eqtyp == "qtelliptic")
            {
              basis = new QTEllipticBasis<2> (order, coeffA, coeffB, coeffC);
            }
          else if (eqtyp == "fowave" || eqtyp == "foqtwave")
            {
              basismats.SetSize (D);
              for (int d = 0; d < D; d++)
                basismats[d] = FOTWaveBasis<2>::Basis (order, d);
              basis = new FOQTWaveBasis<2> (order, coeffA, coeffB);
            }
          else if (eqtyp == "heat")
            {
              basismat = THeatBasis<2>::Basis (order, 0, 0);
            }
          else if (eqtyp == "qtheat")
            {
              basis = new QTHeatBasis<2> (order, coeffA);
            }
          else if (eqtyp == "wave" || eqtyp == "qtwave"
                   || eqtyp == "fowave_reduced")
            {
              basismat = TWaveBasis<2>::Basis (order, basistype,
                                               eqtyp == "fowave_reduced");
              basis = new QTWaveBasis<2> (order, coeffA, coeffB);
            }
          else if (eqtyp == "helmholtz" || eqtyp == "helmholtzconj")
            {
            }
          else
            throw Exception ("TrefftzFESpace::UpdateBasis: unknown eqtyp");
          break;
        }
      case 3:
        {
          if (eqtyp == "laplace")
            basismat = TLapBasis<3>::Basis (order, basistype);
          else if (eqtyp == "qtelliptic")
            {
              basis = new QTEllipticBasis<3> (order, coeffA, coeffB, coeffC);
            }
          else if (eqtyp == "fowave" || eqtyp == "foqtwave")
            {
              basismats.SetSize (D);
              for (int d = 0; d < D; d++)
                basismats[d] = FOTWaveBasis<3>::Basis (order, d);
              basis = new FOQTWaveBasis<3> (order, coeffA, coeffB);
            }
          else if (eqtyp == "heat")
            {
              basismat = THeatBasis<3>::Basis (order, 0, 0);
            }
          else if (eqtyp == "qtheat")
            {
              basis = new QTHeatBasis<3> (order, coeffA);
            }
          else if (eqtyp == "wave" || eqtyp == "qtwave"
                   || eqtyp == "fowave_reduced")
            {
              basismat = TWaveBasis<3>::Basis (order, basistype,
                                               eqtyp == "fowave_reduced");
              basis = new QTWaveBasis<3> (order, coeffA, coeffB);
            }
          else
            throw Exception ("TrefftzFESpace::UpdateBasis: unknown eqtyp");
          break;
        }
      }
  }

  void TrefftzFESpace ::SetCoeff (double acoeff_const)
  {
    coeff_const = acoeff_const;
    UpdateBasis ();
  }

  void TrefftzFESpace ::SetCoeff (shared_ptr<CoefficientFunction> acoeffA,
                                  shared_ptr<CoefficientFunction> acoeffB,
                                  shared_ptr<CoefficientFunction> acoeffC)
  {
    this->coeffA = acoeffA;
    this->coeffB = acoeffB;
    this->coeffC = acoeffC;
    // if (eqtyp.find ("qt") != std::string::npos)
    UpdateBasis ();
  }

  shared_ptr<GridFunction>
  TrefftzFESpace ::GetEWSolution (shared_ptr<CoefficientFunction> acoeffF)
  {
    static Timer t ("QTEll - GetEWSolution");
    RegionTimer reg (t);
    LocalHeap lh (1000 * 1000 * 1000);

    Flags flags;
    flags.SetFlag ("order", order);
    flags.SetFlag ("usescale", this->usescale);
    if (eqtyp == "qtheat" && usescale != 0)
      flags.SetFlag ("usescale", 2);
    shared_ptr<FESpace> fes = make_shared<MonomialFESpace> (ma, flags);
    auto pws = CreateGridFunction (fes, "pws", flags);
    pws->Update ();
    // pws->ConnectAutoUpdate();

    basis->SetRHS (acoeffF);

    // for (auto ei : ma->Elements (VOL))
    ma->IterateElements (VOL, lh, [&] (auto ei, LocalHeap &mlh) {
      HeapReset hr (mlh);
      Array<DofId> dofs;
      fes->GetDofNrs (ei, dofs, VISIBLE_DOF);
      bool hasregdof = false;
      for (DofId d : dofs)
        if (IsRegularDof (d))
          hasregdof = true;
      // assumption here: Either all or no dof is regular
      if (!hasregdof)
        return; // continue;

      FlatVector<double> elvec (dofs.Size (), mlh);

      switch (D)
        {
        case 2:
          {
            Vec<2> scale = 1.0;
            if (eqtyp == "qtheat")
              scale
                  = { ElSize<2> (ei, { 1.0, 0 }), ElSize<2> (ei, { 0, 1.0 }) };
            else if (usescale != 0)
              scale = ElSize<2> (ei);
            basis->GetParticularSolution (ElCenter<2> (ei), scale, elvec, mlh);
            break;
          }
        case 3:
          {
            Vec<3> scale = 1.0;
            if (eqtyp == "qtheat")
              scale = { ElSize<3> (ei, { 1.0, 1.0, 0 }), 0,
                        ElSize<3> (ei, { 0, 0, 1.0 }) };
            else if (usescale != 0)
              scale = ElSize<3> (ei);
            basis->GetParticularSolution (ElCenter<3> (ei), scale, elvec, mlh);
            break;
          }
        }

      pws->SetElementVector (dofs, elvec);
    });

    return pws;
  }

  void TrefftzFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize (0);
    if (!DefinedOn (ei) || ei.VB () != VOL)
      return;
    for (size_t j = ei.Nr () * local_ndof; j < local_ndof * (ei.Nr () + 1);
         j++)
      dnums.Append (j);
  }

  FiniteElement &TrefftzFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement (ei);
    ELEMENT_TYPE eltype = ngel.GetType ();

    if (ei.IsVolume ())
      {
        if (!DefinedOn (ngel))
          {
            return SwitchET (eltype, [&alloc] (auto et) -> FiniteElement & {
              return *new (alloc) ScalarDummyFE<et.ElementType ()> ();
            });
          }

        switch (ma->GetElType (ei))
          {
          case ET_POINT:
          case ET_SEGM:
            {
              // throw Exception ("illegal dim for space-time element");
              break;
            }
          case ET_QUAD:
          case ET_TRIG:
            {
              if (eqtyp == "qtwave")
                {
                  CSR basismat = static_cast<QTWaveBasis<2> *> (basis)->Basis (
                      order, ElCenter<2> (ei));
                  return *(new (alloc) ScalarMappedElement<2> (
                      local_ndof, order, basismat, eltype, ElCenter<2> (ei)));
                }
              else if (eqtyp == "qtelliptic")
                {
                  double scale = 1.0 / ElSize<2> (ei);
                  CSR basismat
                      = static_cast<QTEllipticBasis<2> *> (basis)->Basis (
                          ElCenter<2> (ei), ElSize<2> (ei));
                  return *(new (alloc) ScalarMappedElement<2> (
                      local_ndof, order, basismat, eltype, ElCenter<2> (ei),
                      scale));
                }
              else if (eqtyp == "foqtwave")
                {
                  Vec<2, CSR> qtbasis;
                  for (int d = 0; d < D; d++)
                    qtbasis[d]
                        = static_cast<FOQTWaveBasis<2> *> (basis)->Basis (
                            order, d, ElCenter<2> (ei));
                  return *(new (alloc) BlockMappedElement<2> (
                      local_ndof, order, qtbasis, eltype, ElCenter<2> (ei)));
                }
              else if (eqtyp == "fowave")
                {
                  Vec<2> scale = 1.0;
                  scale[1] = coeff_const;
                  scale = 1.0 / ElSize<2> (ei, scale);
                  scale[1] *= coeff_const;
                  return *(new (alloc) BlockMappedElement<2> (
                      local_ndof, order, basismats, eltype, ElCenter<2> (ei),
                      scale));
                }
              else if (eqtyp == "helmholtz" || eqtyp == "helmholtzconj")
                {
                  return *(new (alloc) PlaneWaveElement<2> (
                      local_ndof, order, eltype, ElCenter<2> (ei),
                      ElSize<2> (ei), coeff_const,
                      (eqtyp == "helmholtz" ? 1 : -1)));
                }
              else if (eqtyp == "heat")
                {
                  Vec<2> scale = 1.0 / sqrt (ElSize<2> (ei));
                  scale[1] = coeff_const * scale[1] * scale[1];
                  return *(new (alloc) ScalarMappedElement<2> (
                      local_ndof, order, basismat, eltype, ElCenter<2> (ei),
                      scale));
                }
              else if (eqtyp == "qtheat")
                {
                  double hx = ElSize<2> (ei, { 1.0, 0 });
                  double ht = ElSize<2> (ei, { 0, 1.0 });
                  Vec<2> scale ({ 1.0 / hx, 1.0 / ht });
                  CSR basismat = static_cast<QTHeatBasis<2> *> (basis)->Basis (
                      ElCenter<2> (ei), hx, ht);
                  return *(new (alloc) ScalarMappedElement<2> (
                      local_ndof, order, basismat, eltype, ElCenter<2> (ei),
                      scale));
                }
              else
                {
                  Vec<2> scale = 1.0;
                  scale[1] = coeff_const;
                  scale = 1.0 / ElSize<2> (ei, scale);
                  scale[1] *= coeff_const;
                  return *(new (alloc) ScalarMappedElement<2> (
                      local_ndof, order, basismat, eltype, ElCenter<2> (ei),
                      scale));
                }
              break;
            }
          case ET_HEX:
          case ET_PRISM:
          case ET_PYRAMID:
          case ET_TET:
            {
              if (eqtyp == "qtwave")
                {
                  CSR basismat = static_cast<QTWaveBasis<3> *> (basis)->Basis (
                      order, ElCenter<3> (ei));
                  return *(new (alloc) ScalarMappedElement<3> (
                      local_ndof, order, basismat, eltype, ElCenter<3> (ei)));
                }
              else if (eqtyp == "qtelliptic")
                {
                  double scale = 1.0 / ElSize<3> (ei);
                  CSR basismat
                      = static_cast<QTEllipticBasis<3> *> (basis)->Basis (
                          ElCenter<3> (ei), ElSize<3> (ei));
                  return *(new (alloc) ScalarMappedElement<3> (
                      local_ndof, order, basismat, eltype, ElCenter<3> (ei),
                      scale));
                }
              else if (eqtyp == "foqtwave")
                {
                  Vec<3, CSR> qtbasis;
                  for (int d = 0; d < D; d++)
                    qtbasis[d]
                        = static_cast<FOQTWaveBasis<3> *> (basis)->Basis (
                            order, d, ElCenter<3> (ei));
                  return *(new (alloc) BlockMappedElement<3> (
                      local_ndof, order, qtbasis, eltype, ElCenter<3> (ei)));
                }
              else if (eqtyp == "fowave")
                {
                  Vec<3> scale = 1.0;
                  scale[2] = coeff_const;
                  scale = 1.0 / ElSize<3> (ei, scale);
                  scale[2] *= coeff_const;
                  return *(new (alloc) BlockMappedElement<3> (
                      local_ndof, order, basismats, eltype, ElCenter<3> (ei),
                      scale));
                }
              else if (eqtyp == "heat")
                {
                  Vec<3> scale = 1.0 / sqrt (ElSize<3> (ei));
                  scale[2] = coeff_const * scale[2] * scale[2];
                  return *(new (alloc) ScalarMappedElement<3> (
                      local_ndof, order, basismat, eltype, ElCenter<3> (ei),
                      scale));
                }
              else if (eqtyp == "qtheat")
                {
                  double hx = ElSize<3> (ei, { 1.0, 1.0, 0 });
                  double ht = ElSize<3> (ei, { 0, 0, 1.0 });
                  Vec<3> scale ({ 1.0 / hx, 1.0 / hx, 1.0 / ht });
                  CSR basismat = static_cast<QTHeatBasis<3> *> (basis)->Basis (
                      ElCenter<3> (ei), hx, ht);
                  return *(new (alloc) ScalarMappedElement<3> (
                      local_ndof, order, basismat, eltype, ElCenter<3> (ei),
                      scale));
                }
              else
                {
                  Vec<3> scale = 1.0;
                  scale[2] = coeff_const;
                  scale = 1.0 / ElSize<3> (ei, scale);
                  scale[2] *= coeff_const;
                  return *(new (alloc) ScalarMappedElement<3> (
                      local_ndof, order, basismat, eltype, ElCenter<3> (ei),
                      scale));
                }
            }
            break;
          }
      }
    // else
    try
      {
        return SwitchET<ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD> (
            eltype, [&alloc] (auto et) -> FiniteElement & {
              return *new (alloc) DummyFE<et.ElementType ()>;
            });
      }
    catch (Exception &e)
      {
        throw Exception ("illegal element type in Trefftz::GetFE");
      }
  }

  DocInfo TrefftzFESpace ::GetDocu ()
  {
    // auto docu = FESpace::GetDocu();
    DocInfo docu;
    docu.short_docu = "Trefftz space for different PDEs. Use kwarg 'eq' to "
                      "choose the PDE, currently implemented are:\n"
                      " - laplace - for Laplace equation\n"
                      " - qtelliptic - for the quasi-Trefftz space for an "
                      "elliptic problem\n"
                      " - wave - for the second order acoustic wave equation\n"
                      " - qtwave - for the quasi-Trefftz space\n"
                      " - fowave - for the first order acoustic wave "
                      "equation, returns TnT (sigv,tauw)\n"
                      " - foqtwave - for the quasi-Trefftz space \n"
                      " - helmholtz - planewaves for the helmholtz equation\n"
                      " - helmholtzconj - returns the complex conjungate of "
                      "the planewaves \n";
    docu.Arg ("eq") = "string\n"
                      "  Choose type of Trefftz functions.";
    docu.Arg ("order") = "int = 1\n"
                         "  Order of finite element space";
    docu.Arg ("dgjumps") = "bool = True\n"
                           "  Enable discontinuous space for DG methods, this "
                           "flag is always True for trefftzfespace.";
    docu.Arg ("complex") = "bool = False\n"
                           "  Set if FESpace should be complex";
    docu.Arg ("useshift") = "bool = True\n"
                            "  shift of basis functins to element center";
    docu.Arg ("usescale") = "bool = True\n"
                            "  scale element basis functions with diam";
    // docu.Arg("useshift") = "bool = True\n"
    //"  use shift of basis functins to element center and scale them";
    // docu.Arg("basistype")
    // docu.Arg("wavespeed")
    return docu;
  }

  static RegisterFESpace<TrefftzFESpace> initi_trefftz ("trefftzfespace");

  //////////////////////////// Trefftz basis ////////////////////////////

  template <int D, typename TFunc>
  void TraversePol (int order, const TFunc &func)
  {
    // traverse polynomials increasing smaller dimensions first
    switch (D)
      {
      case 0:
        break;
      case 1:
        for (int i = 0, ii = 0; i <= order; i++)
          func (ii++, Vec<1, int>{ i });
        break;
      case 2:
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            func (ii++, Vec<2, int>{ j, i });
        break;
      case 3:
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              func (ii++, Vec<3, int>{ k, j, i });
        break;
      case 4:
        for (int i = 0, ii = 0; i <= order; i++)
          for (int j = 0; j <= order - i; j++)
            for (int k = 0; k <= order - i - j; k++)
              for (int l = 0; l <= order - i - j - k; l++)
                func (ii++, Vec<4, int>{ l, k, j, i });
        break;
      default:
        throw Exception ("TraverseDimensions: too many dimensions!");
      }
  }

  template <int D, typename TFunc>
  void TraversePol (Vec<D, size_t> order, const TFunc &func)
  {
    // traverse polynomials increasing smaller dimensions first, up to order in
    // each dimension
    switch (D)
      {
      case 0:
        break;
      case 1:
        for (int i = 0, ii = 0; i <= order[0]; i++)
          func (ii++, Vec<1, int>{ i });
        break;
      case 2:
        for (int i = 0, ii = 0; i <= order[1]; i++)
          for (int j = 0; j <= order[0]; j++)
            func (ii++, Vec<2, int>{ j, i });
        break;
      case 3:
        for (int i = 0, ii = 0; i <= order[2]; i++)
          for (int j = 0; j <= order[1]; j++)
            for (int k = 0; k <= order[0]; k++)
              func (ii++, Vec<3, int>{ k, j, i });
        break;
      case 4:
        for (int i = 0, ii = 0; i <= order[3]; i++)
          for (int j = 0; j <= order[2]; j++)
            for (int k = 0; k <= order[1]; k++)
              for (int l = 0; l <= order[0]; l++)
                func (ii++, Vec<4, int>{ l, k, j, i });
        break;
      default:
        throw Exception ("TraverseDimensions: too many dimensions!");
      }
  }

  template <int D, typename TFunc>
  void TraversePol2 (int order, const TFunc &func)
  {
    // traverse polynomials in order of increasing degree
    switch (D)
      {
      case 0:
        break;
      case 1:
        for (int i = 0, ii = 0; i <= order; i++)
          func (ii++, Vec<1, int>{ i });
        break;
      case 2:
        for (int ord = 0, ii = 0; ord <= order; ord++)
          for (int i = 0; i <= ord; i++)
            {
              int j = ord - i;
              func (ii++, Vec<2, int>{ j, i });
            }
        break;
      case 3:
        for (int ord = 0, ii = 0; ord <= order; ord++)
          for (int i = 0; i <= ord; i++)
            for (int j = 0; j <= ord - i; j++)
              {
                int k = ord - i - j;
                func (ii++, Vec<3, int>{ k, j, i });
              }
        break;
      default:
        throw Exception ("TraverseDimensions: too many dimensions!");
      }
  }

  // k-th coeff of Legendre polynomial of degree n in monomial basis
  constexpr double LegCoeffMonBasis (int n, int k)
  {
    if (n == 0)
      return 1;
    if (k > n)
      return 0;
    if ((n + k) % 2)
      return 0;
    double coeff = pow (2, -n) * pow (-1, floor ((n - k) / 2))
                   * BinCoeff (n, floor ((n - k) / 2)) * BinCoeff (n + k, n);
    // double coeff = pow(2,-n) * pow(-1,k) * BinCoeff(n,k) *
    // BinCoeff(2*n-2*k,n);
    return coeff;
  }

  // k-th coeff of Chebyshev polynomial of degree n in monomial basis
  constexpr double ChebCoeffMonBasis (int n, int k)
  {
    if (n == 0)
      return 1;
    if (k > n)
      return 0;
    if ((n + k) % 2)
      return 0;
    double coeff = pow (2, k - 1) * n * pow (-1, floor ((n - k) / 2))
                   * tgamma ((n + k) / 2)
                   / (tgamma (floor ((n - k) / 2) + 1) * tgamma (k + 1));
    return coeff;
  }

  constexpr int factorial (int n) { return n > 1 ? n * factorial (n - 1) : 1; }

  template <int D> int factorial (Vec<D, int> v)
  {
    int fac = 1;
    for (int i = 0; i < D; i++)
      fac *= factorial (v[i]);
    return fac;
  }

  template <int D>
  CSR TWaveBasis<D>::Basis (int ord, int basistype, int fowave)
  {
    CSR tb;
    const int ndof
        = (BinCoeff (D - 1 + ord, ord) + BinCoeff (D + ord - 2, ord - 1));
    const int npoly = (BinCoeff (D + ord, ord));
    Matrix<> trefftzbasis (ndof, npoly);
    trefftzbasis = 0;
    for (int basis = 0; basis < ndof; basis++)
      {
        int tracker = 0;
        TraversePol<D> (ord, [&] (int i, Vec<D, int> coeff) {
          if (tracker >= 0)
            tracker++;
          int indexmap = PolBasis::IndexMap2<D> (coeff, ord);
          int k = coeff (D - 1);
          if (k == 0 || k == 1)
            {
              switch (basistype)
                {
                case 0:
                  if (tracker > basis)
                    {
                      // trefftzbasis( i, setbasis++ ) = 1.0; //set the l-th
                      // coeff to 1
                      trefftzbasis (basis, indexmap) = 1;
                      tracker = -1;
                    }
                  // i += ndof-1;	//jump to time = 2 if i=0
                  break;
                case 1:
                  if ((k == 0 && basis < BinCoeff (D - 1 + ord, ord))
                      || (k == 1 && basis >= BinCoeff (D - 1 + ord, ord)))
                    {
                      trefftzbasis (basis, indexmap) = 1;
                      for (int exponent : coeff.Range (0, D - 1))
                        trefftzbasis (basis, indexmap)
                            *= LegCoeffMonBasis (basis, exponent);
                    }
                  break;
                case 2:
                  if ((k == 0 && basis < BinCoeff (D - 1 + ord, ord))
                      || (k == 1 && basis >= BinCoeff (D - 1 + ord, ord)))
                    {
                      trefftzbasis (basis, indexmap) = 1;
                      for (int exponent : coeff.Range (0, D - 1))
                        trefftzbasis (basis, indexmap)
                            *= ChebCoeffMonBasis (basis, exponent);
                    }
                  break;
                }
            }
          else if (coeff (D - 1) > 1)
            {
              for (int m = 0; m < D - 1; m++) // rekursive sum
                {
                  Vec<D, int> get_coeff = coeff;
                  get_coeff[D - 1] = get_coeff[D - 1] - 2;
                  get_coeff[m] = get_coeff[m] + 2;
                  trefftzbasis (basis, indexmap)
                      += (coeff (m) + 1) * (coeff (m) + 2)
                         * trefftzbasis (
                             basis, PolBasis::IndexMap2<D> (get_coeff, ord));
                }
              double wavespeed = 1.0;
              trefftzbasis (basis, indexmap)
                  *= wavespeed * wavespeed / (k * (k - 1));
            }
        });
      }
    MatToCSR (trefftzbasis.Rows (fowave, ndof), tb);
    return tb;
  }

  template class TWaveBasis<2>;
  template class TWaveBasis<3>;
  template class TWaveBasis<4>;

  template <int D>
  CSR THeatBasis<D>::Basis (int ord, int basistype, int fowave)
  {
    CSR tb;
    const int ndof = BinCoeff (D - 1 + ord, ord);
    const int npoly = BinCoeff (D + ord, ord);
    Matrix<> trefftzbasis (ndof, npoly);
    trefftzbasis = 0;
    for (int basis = 0; basis < ndof; basis++)
      {
        int tracker = 0;
        TraversePol<D> (ord, [&] (int i, Vec<D, int> coeff) {
          int indexmap = PolBasis::IndexMap2<D> (coeff, ord);
          int sum = 0;
          for (int m = 0; m <= D - 1; m++)
            sum += coeff[m];
          if (coeff[D - 1] == 0 && tracker >= basis)
            {
              trefftzbasis (basis, indexmap) = 1;
              tracker = -1;
            }
          else if (coeff[D - 1] == 0 && tracker >= 0)
            {
              tracker++;
            }
          else if (coeff[D - 1] > 0 && coeff[D - 1] <= ord / 2.0 && sum < ord)
            {
              for (int m = 0; m < D - 1; m++) // rekursive sum
                {
                  Vec<D, int> get_coeff = coeff;
                  get_coeff[D - 1] = get_coeff[D - 1] - 1;
                  get_coeff[m] = get_coeff[m] + 2;
                  trefftzbasis (basis, indexmap)
                      += (coeff[m] + 1) * (coeff[m] + 2)
                         * trefftzbasis (
                             basis, PolBasis::IndexMap2<D> (get_coeff, ord));
                }
              trefftzbasis (basis, indexmap) *= 1.0 / coeff[D - 1];
            }
        });
      }
    MatToCSR (trefftzbasis.Rows (fowave, ndof), tb);
    return tb;
  }

  template class THeatBasis<2>;
  template class THeatBasis<3>;

  template <int D> CSR TLapBasis<D>::Basis (int ord, int basistype)
  {
    CSR tb;
    const int ndof
        = (BinCoeff (D - 1 + ord, ord) + BinCoeff (D - 1 + ord - 1, ord - 1));
    const int npoly = (BinCoeff (D + ord, ord));
    Matrix<> trefftzbasis (ndof, npoly);
    trefftzbasis = 0;
    for (int basis = 0; basis < ndof; basis++)
      {
        int tracker = 0;
        TraversePol<D> (ord, [&] (int i, Vec<D, int> coeff) {
          if (tracker >= 0)
            tracker++;
          int indexmap = PolBasis::IndexMap2<D> (coeff, ord);
          int k = coeff (D - 1);
          if (k == 0 || k == 1)
            {
              if (tracker > basis)
                {
                  // trefftzbasis( i, setbasis++ ) = 1.0; //set the l-th coeff
                  // to 1
                  trefftzbasis (basis, indexmap) = 1;
                  tracker = -1;
                }
            }
          else if (coeff (D - 1) > 1)
            {
              for (int m = 0; m < D - 1; m++) // rekursive sum
                {
                  Vec<D, int> get_coeff = coeff;
                  get_coeff[D - 1] = get_coeff[D - 1] - 2;
                  get_coeff[m] = get_coeff[m] + 2;
                  trefftzbasis (basis, indexmap)
                      -= (coeff (m) + 1) * (coeff (m) + 2)
                         * trefftzbasis (
                             basis, PolBasis::IndexMap2<D> (get_coeff, ord));
                }
              double lapcoeff = 1.0;
              trefftzbasis (basis, indexmap)
                  *= lapcoeff * lapcoeff / (k * (k - 1));
            }
        });
      }
    MatToCSR (trefftzbasis, tb);
    return tb;
  }

  template class TLapBasis<1>;
  template class TLapBasis<2>;
  template class TLapBasis<3>;

  template <int D> CSR FOTWaveBasis<D>::Basis (int ord, int rdim)
  {
    const int ndof = D * BinCoeff (ord + D - 1, D - 1);
    const int npoly = BinCoeff (D + ord, ord);
    Array<Matrix<>> trefftzbasis (D);
    for (int d = 0; d < D; d++)
      {
        trefftzbasis[d].SetSize (ndof, npoly);
        trefftzbasis[d] = 0;
      }
    for (int basis = 0; basis < ndof; basis++)
      {
        int tracker = 0;
        TraversePol<D> (ord, [&] (int i, Vec<D, int> coeff) {
          if (tracker >= 0)
            tracker++;
          int indexmap = PolBasis::IndexMap2<D> (coeff, ord);
          if (coeff (D - 1) == 0 && tracker * (D) > basis)
            {
              int d = basis % (D);
              trefftzbasis[d](basis, indexmap) = 1;
              tracker = -1;
            }
          else if (coeff (D - 1) > 0)
            {
              int k = coeff (D - 1);
              for (int d = 0; d < D - 1; d++)
                {
                  Vec<D, int> get_coeff = coeff;
                  get_coeff[D - 1] = get_coeff[D - 1] - 1;
                  get_coeff[d] = get_coeff[d] + 1;
                  trefftzbasis[d](basis, indexmap)
                      = (-1.0 / k) * (coeff (d) + 1)
                        * trefftzbasis[D - 1](
                            basis, PolBasis::IndexMap2<D> (get_coeff, ord));
                  trefftzbasis[D - 1](basis, indexmap)
                      += (-1.0 / k) * (coeff (d) + 1)
                         * trefftzbasis[d](
                             basis, PolBasis::IndexMap2<D> (get_coeff, ord));
                }
            }
        });
      }
    Array<CSR> tb (D);
    for (int d = 0; d < D; d++)
      {
        // cout << d << endl << trefftzbasis[d] << endl;
        MatToCSR (trefftzbasis[d], tb[d]);
      }
    return tb[rdim];
  }

  template class FOTWaveBasis<2>;
  template class FOTWaveBasis<3>;
  template class FOTWaveBasis<4>;

  //////////////////////////// quasi-Trefftz basis ////////////////////////////

  template <int D>
  CSR QTEllipticBasis<D>::Basis (Vec<D> ElCenter, double elsize)
  {
    lock_guard<mutex> lock (gentrefftzbasis);
    int order = this->order;
    string encode = to_string (order) + to_string (elsize);
    for (int i = 0; i < D; i++)
      encode += to_string (ElCenter[i]);

    if (gtbstore[encode][0].Size () == 0)
      {
        IntegrationPoint ip (ElCenter, 0);
        Mat<D, D> dummy;
        FE_ElementTransformation<D, D> et (D == 3   ? ET_TET
                                           : D == 2 ? ET_TRIG
                                                    : ET_SEGM,
                                           dummy);
        MappedIntegrationPoint<D, D> mip (ip, et, 0);
        for (int i = 0; i < D; i++)
          mip.Point ()[i] = ElCenter[i];

        const int ndiffs = (BinCoeff (D + order - 1, order - 1));
        Vector<Matrix<>> AA (ndiffs);
        Vector<Vector<>> BB (ndiffs);
        Vector<> CC (ndiffs);

        TraversePol<D> (order - 1, [&] (int i, Vec<D, int> coeff) {
          int index = PolBasis::IndexMap2<D> (coeff, order - 1);
          AA[index].SetSize (D, D);
          BB[index].SetSize (D);
          AAder[index]->Evaluate (mip, AA[index].AsVector ());
          BBder[index]->Evaluate (mip, BB[index]);
          CC[index] = CCder[index]->Evaluate (mip);
        });

        const int ndof = (BinCoeff (D - 1 + order, order)
                          + BinCoeff (D - 1 + order - 1, order - 1));
        const int npoly = (BinCoeff (D + order, order));
        Matrix<> qtbasis (ndof, npoly);
        qtbasis = 0;
        // init qtbasis
        TraversePol<D> (order, [&] (int i, Vec<D, int> coeff) {
          if (coeff[D - 1] > 1)
            return;
          int indexmap = PolBasis::IndexMap2<D> (coeff, order);
          qtbasis (i, indexmap) = 1;
        });
        // start recursion
        TraversePol2<D> (order, [&] (int i, Vec<D, int> coeff) {
          if (coeff (D - 1) <= 1)
            return;
          int indexmap = PolBasis::IndexMap2<D> (coeff, order);
          Vec<D, int> mii = coeff;
          mii[D - 1] = mii[D - 1] - 2;
          for (int j = 0; j < D; j++)
            {
              Vec<D, int> ej = 0;
              ej[j] = 1;

              TraversePol<D> (mii + ej, [&] (int i2, Vec<D, int> mil) {
                // matrix coeff A
                for (int m = 0; m < D; m++)
                  {
                    if (i2 == 0 && m == D - 1 && j == D - 1)
                      continue;
                    Vec<D, int> em = 0;
                    em[m] = 1;
                    qtbasis.Col (indexmap)
                        -= factorial (mii + ej) / factorial (mil)
                           * (AA[IndexMap2<D> (mil, order - 1)]) (j, m)
                           * pow (elsize, vsum<D, int> (mil))
                           * (mii[m] + ej[m] - mil[m] + 1)
                           * qtbasis.Col (PolBasis::IndexMap2<D> (
                               mii + ej - mil + em, order));
                  }
                // vec coeff B
                qtbasis.Col (indexmap)
                    += factorial (mii + ej) / factorial (mil)
                       * (BB[IndexMap2<D> (mil, order - 1)]) (j)*pow (
                           elsize, vsum<D, int> (mil) + 1)
                       * qtbasis.Col (
                           PolBasis::IndexMap2<D> (mii + ej - mil, order));

                // scal coeff C
                if (j == 0 && mil[0] <= mii[0])
                  qtbasis.Col (indexmap)
                      += factorial (mii) / factorial (mil)
                         * CC[IndexMap2<D> (mil, order - 1)]
                         * pow (elsize, vsum<D, int> (mil) + 2)
                         * qtbasis.Col (
                             PolBasis::IndexMap2<D> (mii - mil, order));
              });
            }
          Vec<D, int> eD = 0;
          eD[D - 1] = 2;
          qtbasis.Col (indexmap)
              *= 1.0 / factorial (mii + eD) / (AA[0](D - 1, D - 1));
        });
        MatToCSR (qtbasis, gtbstore[encode]);
      }

    if (gtbstore[encode].Size () == 0)
      {
        stringstream str;
        str << "failed to generate trefftz basis of order " << order << endl;
        throw Exception (str.str ());
      }

    return gtbstore[encode];
    // CSR tb;
    // MatToCSR (qtbasis, tb);
    // return tb;
  }

  template <int D>
  void
  QTEllipticBasis<D>::GetParticularSolution (Vec<D> ElCenter, Vec<D> elsize,
                                             FlatVector<> sol, LocalHeap &lh)
  {
    double hx = elsize[0];
    static Timer t ("QTEll - GetParticularSolution");
    RegionTimer reg (t);
    IntegrationPoint ip (ElCenter, 0);
    Mat<D, D> dummy;
    FE_ElementTransformation<D, D> et (D == 3   ? ET_TET
                                       : D == 2 ? ET_TRIG
                                                : ET_SEGM,
                                       dummy);
    MappedIntegrationPoint<D, D> mip (ip, et, 0);
    for (int i = 0; i < D; i++)
      mip.Point ()[i] = ElCenter[i];

    int ndiffs = (BinCoeff (D + order - 1, order - 1));
    FlatVector<Matrix<>> AA (ndiffs, lh);
    FlatVector<Vector<>> BB (ndiffs, lh);
    FlatVector<> CC (ndiffs, lh);
    ndiffs = (BinCoeff (D + order, order));
    FlatVector<> FF (ndiffs, lh);

    TraversePol<D> (order, [&] (int i, Vec<D, int> coeff) {
      int index = PolBasis::IndexMap2<D> (coeff, order);
      FF[index] = FFder[index]->Evaluate (mip);
      if (vsum<D, int> (coeff) < order)
        {
          index = PolBasis::IndexMap2<D> (coeff, order - 1);
          AA[index].AssignMemory (D, D, lh);
          BB[index].AssignMemory (D, lh);
          AAder[index]->Evaluate (mip, AA[index].AsVector ());
          BBder[index]->Evaluate (mip, BB[index]);
          CC[index] = CCder[index]->Evaluate (mip);
        }
    });

    sol = 0;
    // start recursion
    TraversePol2<D> (order, [&] (int i, Vec<D, int> mii) {
      if (mii (D - 1) <= 1)
        return;
      int indexmap = PolBasis::IndexMap2<D> (mii, order);
      mii[D - 1] = mii[D - 1] - 2;
      for (int j = 0; j < D; j++)
        {
          Vec<D, int> ej = 0;
          ej[j] = 1;

          TraversePol<D> (mii + ej, [&] (int i2, Vec<D, int> mil) {
            // matrix coeff A
            for (int m = 0; m < D; m++)
              {
                if (i2 == 0 && m == D - 1 && j == D - 1)
                  continue;
                Vec<D, int> em = 0;
                em[m] = 1;
                sol (indexmap) -= factorial (mii + ej) / factorial (mil)
                                  * (AA[IndexMap2<D> (mil, order - 1)]) (j, m)
                                  * pow (hx, vsum<D, int> (mil))
                                  * (mii[m] + ej[m] - mil[m] + 1)
                                  * sol (PolBasis::IndexMap2<D> (
                                      mii + ej - mil + em, order));
              }
            // vec coeff B
            sol (indexmap)
                += factorial (mii + ej) / factorial (mil)
                   * (BB[IndexMap2<D> (mil, order - 1)]) (j)*pow (
                       hx, vsum<D, int> (mil) + 1)
                   * sol (PolBasis::IndexMap2<D> (mii + ej - mil, order));

            // scal coeff C
            if (j == 0 && mil[0] <= mii[0])
              sol (indexmap)
                  += factorial (mii) / factorial (mil)
                     * CC[IndexMap2<D> (mil, order - 1)]
                     * pow (hx, vsum<D, int> (mil) + 2)
                     * sol (PolBasis::IndexMap2<D> (mii - mil, order));
          });
        }
      sol (indexmap) += -FF[PolBasis::IndexMap2<D> (mii, order)]
                        * pow (hx, vsum<D, int> (mii) + 2);
      Vec<D, int> eD = 0;
      eD[D - 1] = 2;
      sol (indexmap) *= 1.0 / factorial (mii + eD) / (AA[0](D - 1, D - 1));
    });
  }

  template class QTEllipticBasis<1>;
  template class QTEllipticBasis<2>;
  template class QTEllipticBasis<3>;

  template <int D>
  CSR QTWaveBasis<D>::Basis (int ord, Vec<D> ElCenter, double elsize,
                             int basistype)
  {
    lock_guard<mutex> lock (gentrefftzbasis);
    string encode = to_string (ord) + to_string (elsize);
    for (int i = 0; i < D - 1; i++)
      encode += to_string (ElCenter[i]);

    if (gtbstore[encode][0].Size () == 0)
      {
        IntegrationPoint ip (ElCenter, 0);
        Mat<D - 1, D - 1> dummy;
        FE_ElementTransformation<D - 1, D - 1> et (D == 4   ? ET_TET
                                                   : D == 3 ? ET_TRIG
                                                            : ET_SEGM,
                                                   dummy);
        MappedIntegrationPoint<D - 1, D - 1> mip (ip, et, 0);
        for (int i = 0; i < D - 1; i++)
          mip.Point ()[i] = ElCenter[i];

        Matrix<> BB (ord, (ord - 1) * (D == 3) + 1);
        Matrix<> AA (ord - 1, (ord - 2) * (D == 3) + 1);

        TraversePol<D - 1> (order - 1, [&] (int i, Vec<D - 1, int> coeff) {
          int nx = coeff[0];
          int ny = D > 2 ? coeff[1] : 0;
          double fac = (factorial (nx) * factorial (ny));
          int index = PolBasis::IndexMap2<D - 1> (coeff, order - 1);
          BB (nx, ny)
              = BBder[index]->Evaluate (mip) / fac * pow (elsize, nx + ny);
          if (vsum<D - 1, int> (coeff) < ord - 1)
            {
              index = PolBasis::IndexMap2<D - 1> (coeff, order - 2);
              AA (nx, ny)
                  = AAder[index]->Evaluate (mip) / fac * pow (elsize, nx + ny);
            }
        });

        const int ndof
            = (BinCoeff (D - 1 + ord, ord) + BinCoeff (D + ord - 2, ord - 1));
        const int npoly = BinCoeff (D + ord, ord);
        Matrix<> qtbasis (ndof, npoly);
        qtbasis = 0;

        for (int t = 0, basisn = 0; t < 2; t++)
          for (int x = 0; x <= ord - t; x++)
            for (int y = 0; y <= (ord - x - t) * (D == 3); y++)
              {
                Vec<D, int> index;
                index[D - 1] = t;
                index[0] = x;
                if (D == 3)
                  index[1] = y;
                qtbasis (basisn++, PolBasis::IndexMap2<D> (index, ord)) = 1.0;
              }

        for (int basisn = 0; basisn < ndof; basisn++)
          {
            for (int ell = 0; ell < ord - 1; ell++)
              {
                for (int t = 0; t <= ell; t++)
                  {
                    for (int x = (D == 2 ? ell - t : 0); x <= ell - t; x++)
                      {
                        int y = ell - t - x;
                        Vec<D, int> index;
                        index[1] = y;
                        index[0] = x;
                        index[D - 1] = t + 2;
                        double *newcoeff = &qtbasis (
                            basisn, PolBasis::IndexMap2<D> (index, ord));
                        *newcoeff = 0;

                        for (int betax = 0; betax <= x; betax++)
                          for (int betay = (D == 3) ? 0 : y; betay <= y;
                               betay++)
                            {
                              index[1] = betay;
                              index[0] = betax + 1;
                              index[D - 1] = t;
                              int getcoeffx
                                  = PolBasis::IndexMap2<D> (index, ord);
                              index[1] = betay + 1;
                              index[0] = betax;
                              index[D - 1] = t;
                              int getcoeffy
                                  = PolBasis::IndexMap2<D> (index, ord);
                              index[1] = betay;
                              index[0] = betax + 2;
                              index[D - 1] = t;
                              int getcoeffxx
                                  = PolBasis::IndexMap2<D> (index, ord);
                              index[1] = betay + 2;
                              index[0] = betax;
                              index[D - 1] = t;
                              int getcoeffyy
                                  = PolBasis::IndexMap2<D> (index, ord);

                              *newcoeff
                                  += (betax + 2) * (betax + 1)
                                         / ((t + 2) * (t + 1) * AA (0))
                                         * BB (x - betax, y - betay)
                                         * qtbasis (basisn, getcoeffxx)
                                     + (x - betax + 1) * (betax + 1)
                                           / ((t + 2) * (t + 1) * AA (0))
                                           * BB (x - betax + 1, y - betay)
                                           * qtbasis (basisn, getcoeffx);
                              if (D == 3)
                                *newcoeff
                                    += (betay + 2) * (betay + 1)
                                           / ((t + 2) * (t + 1) * AA (0))
                                           * BB (x - betax, y - betay)
                                           * qtbasis (basisn, getcoeffyy)
                                       + (y - betay + 1) * (betay + 1)
                                             / ((t + 2) * (t + 1) * AA (0))
                                             * BB (x - betax, y - betay + 1)
                                             * qtbasis (basisn, getcoeffy);
                              if (betax + betay == x + y)
                                continue;
                              index[1] = betay;
                              index[0] = betax;
                              index[D - 1] = t + 2;
                              int getcoeff
                                  = PolBasis::IndexMap2<D> (index, ord);

                              *newcoeff -= AA (x - betax, y - betay)
                                           * qtbasis (basisn, getcoeff)
                                           / AA (0);
                            }
                      }
                  }
              }
          }

        MatToCSR (qtbasis, gtbstore[encode]);
      }

    if (gtbstore[encode].Size () == 0)
      {
        stringstream str;
        str << "failed to generate trefftz basis of order " << ord << endl;
        throw Exception (str.str ());
      }

    return gtbstore[encode];
  }

  template class QTWaveBasis<2>;
  template class QTWaveBasis<3>;

  template <int D>
  CSR FOQTWaveBasis<D>::Basis (int ord, int rdim, Vec<D> ElCenter,
                               double elsize)
  {
    lock_guard<mutex> lock (gentrefftzbasis);
    string encode = to_string (ord) + to_string (elsize);
    for (int i = 0; i < D - 1; i++)
      encode += to_string (ElCenter[i]);

    if (gtbstore[0][encode][0].Size () == 0)
      {
        IntegrationPoint ip (ElCenter, 0);
        Mat<D - 1, D - 1> dummy;
        FE_ElementTransformation<D - 1, D - 1> et (D == 4   ? ET_TET
                                                   : D == 3 ? ET_TRIG
                                                            : ET_SEGM,
                                                   dummy);
        MappedIntegrationPoint<D - 1, D - 1> mip (ip, et, 0);
        for (int i = 0; i < D - 1; i++)
          mip.Point ()[i] = ElCenter[i];

        Matrix<> BB (ord, (ord - 1) * (D == 3) + 1);
        Matrix<> AA (ord, (ord - 1) * (D == 3) + 1);

        TraversePol<D - 1> (order - 1, [&] (int i, Vec<D - 1, int> coeff) {
          int nx = coeff[0];
          int ny = D > 2 ? coeff[1] : 0;
          double fac = (factorial (nx) * factorial (ny));
          int index = PolBasis::IndexMap2<D - 1> (coeff, order - 1);
          BB (nx, ny)
              = BBder[index]->Evaluate (mip) / fac * pow (elsize, nx + ny);
          AA (nx, ny)
              = AAder[index]->Evaluate (mip) / fac * pow (elsize, nx + ny);
        });

        const int ndof = D * BinCoeff (ord + D - 1, D - 1);
        const int npoly = BinCoeff (D + ord, ord);
        Array<Matrix<>> qtbasis (D);
        for (int d = 0; d < D; d++)
          {
            qtbasis[d].SetSize (ndof, npoly);
            qtbasis[d] = 0;
          }

        for (int d = 0, basisn = 0; d < D; d++)
          {
            for (int x = 0; x <= ord; x++)
              for (int y = 0; y <= (ord - x) * (D == 3); y++)
                {
                  Vec<D, int> index;
                  index[1] = y;
                  index[0] = x;
                  index[D - 1] = 0;
                  qtbasis[d](basisn++, PolBasis::IndexMap2<D> (index, ord))
                      = 1.0;
                }
          }

        for (int basisn = 0; basisn < ndof; basisn++)
          {
            for (int ell = 0; ell < ord; ell++)
              {
                for (int t = 0; t <= ell; t++)
                  {
                    for (int x = (D == 2 ? ell - t : 0); x <= ell - t; x++)
                      {
                        int y = ell - t - x;
                        Vec<D, int> index;
                        index[1] = y;
                        index[0] = x;
                        index[D - 1] = t + 1;
                        int newindex = PolBasis::IndexMap2<D> (index, ord);
                        double *newcoefft = &qtbasis[D - 1](basisn, newindex);
                        for (int d = 0; d < D - 1; d++)
                          {
                            double *newcoeff = &qtbasis[d](basisn, newindex);

                            index[1] = y + (d == 1);
                            index[0] = x + (d == 0);
                            index[D - 1] = t;
                            int getcoeff = PolBasis::IndexMap2<D> (index, ord);
                            *newcoeff = -qtbasis[D - 1](basisn, getcoeff)
                                        * index[d] / (t + 1) / BB (0);
                            *newcoefft -= qtbasis[d](basisn, getcoeff)
                                          * index[d] / (t + 1) / AA (0);
                            for (int betax = 0; betax <= x; betax++)
                              for (int betay = (D == 2) ? y : 0; betay <= y;
                                   betay++)
                                {
                                  if (betax + betay == x + y)
                                    continue;
                                  index[1] = betay;
                                  index[0] = betax;
                                  index[D - 1] = t + 1;
                                  int getcoeff
                                      = PolBasis::IndexMap2<D> (index, ord);
                                  *newcoeff -= BB (x - betax, y - betay)
                                               * qtbasis[d](basisn, getcoeff)
                                               / BB (0);
                                  if (d == 0)
                                    *newcoefft
                                        -= AA (x - betax, y - betay)
                                           * qtbasis[D - 1](basisn, getcoeff)
                                           / AA (0);
                                }
                          }
                      }
                  }
              }
          }

        for (int d = 0; d < D; d++)
          {
            MatToCSR (qtbasis[d], gtbstore[d][encode]);
          }
      }

    if (gtbstore[0][encode].Size () == 0)
      {
        stringstream str;
        str << "failed to generate trefftz basis of order " << ord << endl;
        throw Exception (str.str ());
      }

    return gtbstore[rdim][encode];
  }

  template class FOQTWaveBasis<2>;
  template class FOQTWaveBasis<3>;

  template <int D>
  CSR QTHeatBasis<D>::Basis (Vec<D> ElCenter, double hx, double ht)
  {
    // lock_guard<mutex> lock (gentrefftzbasis);
    // int order = this->order;
    // string encode = to_string (order) + to_string (hx) + to_string (ht);
    // for (int i = 0; i < D-1; i++)
    // encode += to_string (ElCenter[i]);

    // if (gtbstore[encode][0].Size () == 0)
    {
      IntegrationPoint ip (ElCenter, 0);
      Mat<D, D> dummy;
      FE_ElementTransformation<D, D> et (D == 3 ? ET_TET : ET_TRIG, dummy);
      MappedIntegrationPoint<D, D> mip (ip, et, 0);
      for (int i = 0; i < D; i++)
        mip.Point ()[i] = ElCenter[i];

      const int ndiffs = (BinCoeff (D + order - 1, order - 1));
      Vector<Matrix<>> AA (ndiffs);

      TraversePol<D> (order - 1, [&] (int i, Vec<D, int> coeff) {
        int index = PolBasis::IndexMap2<D> (coeff, order - 1);
        AA[index].SetSize (D - 1, D - 1);
        AAder[index]->Evaluate (mip, AA[index].AsVector ());
      });

      const int ndof = (BinCoeff (D - 1 + order, order)
                        + BinCoeff (D - 1 + order - 1, order - 1));
      const int npoly = (BinCoeff (D + order, order));
      Matrix<> qtbasis (ndof, npoly);
      qtbasis = 0;
      // init qtbasis
      int counter = 0;
      TraversePol<D> (order, [&] (int i, Vec<D, int> coeff) {
        if (coeff[0] > 1)
          return;
        int indexmap = PolBasis::IndexMap2<D> (coeff, order);
        qtbasis (counter++, indexmap) = 1;
      });
      // start recursion
      TraversePol2<D> (order, [&] (int i, Vec<D, int> coeff) {
        if (coeff[0] <= 1)
          return;
        int indexmap = PolBasis::IndexMap2<D> (coeff, order);
        Vec<D, int> mii = coeff;
        mii[0] = mii[0] - 2;
        // mii[D - 1] = mii[D - 1] + 1;
        Vec<D, int> et = 0;
        et[D - 1] = 1;
        qtbasis.Col (indexmap)
            += hx * hx / ht * (mii[D - 1] + 1) / coeff[0] / (coeff[0] - 1)
               * qtbasis.Col (PolBasis::IndexMap2<D> (mii + et, order));
        for (int j = 0; j < D - 1; j++)
          {
            Vec<D, int> ej = 0;
            ej[j] = 1;
            TraversePol<D> (mii + ej, [&] (int i2, Vec<D, int> mil) {
              // matrix coeff A
              for (int m = 0; m < D - 1; m++)
                {
                  if (i2 == 0 && m == 0 && j == 0)
                    continue;
                  Vec<D, int> em = 0;
                  em[m] = 1;
                  qtbasis.Col (indexmap)
                      -= (AA[IndexMap2<D> (mil, order - 1)]) (j, m)
                         * pow (hx, vsum<D - 1, int> (mil))
                         * pow (ht, mil[D - 1]) * (mii[m] + ej[m] - mil[m] + 1)
                         * (mii[j] + 1) / (mii[0] + 2) / (mii[0] + 1)
                         / factorial (mil)
                         * qtbasis.Col (PolBasis::IndexMap2<D> (
                             mii + ej - mil + em, order));
                }
            });
          }
        qtbasis.Col (indexmap) *= 1.0 / (AA[0](0, 0));
      });

      // MatToCSR (qtbasis, gtbstore[encode]);
      CSR tb;
      MatToCSR (qtbasis, tb);
      return tb;
    }

    // if (gtbstore[encode].Size () == 0)
    //{
    // stringstream str;
    // str << "failed to generate trefftz basis of order " << order << endl;
    // throw Exception (str.str ());
    //}

    // return gtbstore[encode];
  }

  template <int D>
  void QTHeatBasis<D>::GetParticularSolution (Vec<D> ElCenter, Vec<D> elsize,
                                              FlatVector<> sol, LocalHeap &lh)
  {
    double hx = elsize[0];
    double ht = elsize[D - 1];
    IntegrationPoint ip (ElCenter, 0);
    Mat<D, D> dummy;
    FE_ElementTransformation<D, D> et (D == 3 ? ET_TET : ET_TRIG, dummy);
    MappedIntegrationPoint<D, D> mip (ip, et, 0);
    for (int i = 0; i < D; i++)
      mip.Point ()[i] = ElCenter[i];

    int ndiffs = (BinCoeff (D + order - 1, order - 1));
    Vector<Matrix<>> AA (ndiffs);
    ndiffs = (BinCoeff (D + order, order));
    FlatVector<> FF (ndiffs, lh);

    TraversePol<D> (order, [&] (int i, Vec<D, int> coeff) {
      int index = PolBasis::IndexMap2<D> (coeff, order);
      FF[index] = FFder[index]->Evaluate (mip);
      if (vsum<D, int> (coeff) < order)
        {
          index = PolBasis::IndexMap2<D> (coeff, order - 1);
          AA[index].SetSize (D - 1, D - 1);
          AAder[index]->Evaluate (mip, AA[index].AsVector ());
        }
    });

    const int ndof = (BinCoeff (D - 1 + order, order)
                      + BinCoeff (D - 1 + order - 1, order - 1));
    const int npoly = (BinCoeff (D + order, order));
    sol = 0;
    // start recursion
    TraversePol2<D> (order, [&] (int i, Vec<D, int> coeff) {
      if (coeff[0] <= 1)
        return;
      int indexmap = PolBasis::IndexMap2<D> (coeff, order);
      Vec<D, int> mii = coeff;
      mii[0] = mii[0] - 2;
      Vec<D, int> et = 0;
      et[D - 1] = 1;
      sol (indexmap) += hx * hx / ht * (mii[D - 1] + 1) / coeff[0]
                        / (coeff[0] - 1)
                        * sol (PolBasis::IndexMap2<D> (mii + et, order));
      for (int j = 0; j < D - 1; j++)
        {
          Vec<D, int> ej = 0;
          ej[j] = 1;
          TraversePol<D> (mii + ej, [&] (int i2, Vec<D, int> mil) {
            // matrix coeff A
            for (int m = 0; m < D - 1; m++)
              {
                if (i2 == 0 && m == 0 && j == 0)
                  continue;
                Vec<D, int> em = 0;
                em[m] = 1;
                sol (indexmap)
                    -= (AA[IndexMap2<D> (mil, order - 1)]) (j, m)
                       * pow (hx, vsum<D - 1, int> (mil))
                       * pow (ht, mil[D - 1]) * (mii[m] + ej[m] - mil[m] + 1)
                       * (mii[j] + 1) / (mii[0] + 2) / (mii[0] + 1)
                       / factorial (mil)
                       * sol (PolBasis::IndexMap2<D> (mii + ej - mil + em,
                                                      order));
              }
          });
        }
      sol (indexmap) += -FF[PolBasis::IndexMap2<D> (mii, order)]
                        * pow (hx, vsum<D - 1, int> (coeff))
                        * pow (ht, coeff[D - 1]) / factorial (coeff);
      sol (indexmap) *= 1.0 / (AA[0](0, 0));
    });
    // cout << sol << endl;
  }

  // template class QTHeatBasis<1>;
  template class QTHeatBasis<2>;
  template class QTHeatBasis<3>;
}

#ifdef NGS_PYTHON
void ExportTrefftzFESpace (py::module m)
{
  using namespace ngcomp;

  ExportFESpace<TrefftzFESpace> (m, "trefftzfespace")
      .def ("GetDocu", &TrefftzFESpace::GetDocu)
      .def ("GetNDof", &TrefftzFESpace::GetNDof)
      .def ("SetCoeff",
            static_cast<void (TrefftzFESpace::*) (double)> (
                &TrefftzFESpace::SetCoeff),
            py::arg ("coeff_const"))
      .def (
          "SetCoeff",
          static_cast<void (TrefftzFESpace::*) (
              shared_ptr<CoefficientFunction>, shared_ptr<CoefficientFunction>,
              shared_ptr<CoefficientFunction>)> (&TrefftzFESpace::SetCoeff),
          R"mydelimiter(
                Set coefficient of Trefftz space. 

                For an elliptic problem, the coefficients are given by
                - div(coeffA*grad(u)) + coeffB*grad(u) + coeffC u = 0

                For the first order wave equation, the coefficients are given by
                grad(v) + coeffB dt sigma = 0
                div(sigma) + 1/coeffA**2 dt v = 0

                For the second order wave equation, the coefficients are given by
                - div(1/coeffB grad(u)) + 1/coeffA**2 dtt u = 0

                Parameters
                ----------
                coeffA : CoefficientFunction
                    Coefficient A
                coeffB : CoefficientFunction
                    Coefficient B
                coeffC : CoefficientFunction
                    Coefficient C
            )mydelimiter",
          py::arg ("acoeffA"), py::arg ("acoeffB") = nullptr,
          py::arg ("acoeffC") = nullptr)
      .def ("GetEWSolution",
            static_cast<shared_ptr<GridFunction> (TrefftzFESpace::*) (
                shared_ptr<CoefficientFunction>)> (
                &TrefftzFESpace::GetEWSolution),
            R"mydelimiter(
                Compute a element-wise particular solution for given right hand side.

                Parameters
                ----------
                coeffF : CoefficientFunction
                    Right hand side
            )mydelimiter",
            py::arg ("acoeffF"));

  // ExportFESpace<FOTWaveFESpace, CompoundFESpace> (m, "FOTWave");
}
#endif // NGS_PYTHON
