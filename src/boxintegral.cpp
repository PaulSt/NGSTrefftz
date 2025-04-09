#include "boxintegral.hpp"

#include <python_comp.hpp>

using namespace ngfem;
using namespace ngcomp;
using namespace ngbla;

bool symbolic_integrator_uses_diff = false;

template <int D>
INLINE double GetElementSizeCenter (const ElementTransformation &trafo,
                                    FlatVector<> &center, LocalHeap &lh)
{
  center = 0;
  auto &ip = *(
      new (lh) IntegrationPoint (1.0 / (D + 1), 1.0 / (D + 1), 1.0 / (D + 1)));
  auto &mip = *(new (lh) MappedIntegrationPoint<D, D> (ip, trafo));
  center = mip.GetPoint ();
  if (D == 1)
    return mip.GetJacobiDet ();
  else if (D == 2)
    return sqrt (mip.GetJacobiDet ());
  else
    return cbrt (mip.GetJacobiDet ());
}

/**
 * Maps reference points to physical points and computes weights.
 *
 * @param trafo The element transformation.
 * @param pts_in The input reference points.
 * @param pts_out The output physical points.
 * @param wts_in The input weights.
 * @param wts_out The output weights.
 * @param lh The LocalHeap object.
 */
template <int D>
INLINE void
MapRefPoints (const ElementTransformation &trafo,
              FlatMatrixFixWidth<D> &pts_in, FlatMatrixFixWidth<D> &pts_out,
              FlatVector<> &wts_in, FlatVector<> &wts_out, double box_length,
              bool scale_with_elsize, LocalHeap &lh, VorB element_vb = VOL)
{
  FlatVector<> center (D, lh);
  double h = GetElementSizeCenter<D> (trafo, center, lh);
  if (scale_with_elsize)
    box_length *= h;
  const double box_factor = pow (box_length, element_vb == VOL ? D : D - 1);
  for (size_t i = 0; i < pts_in.Height (); i++)
    {
      for (int d = 0; d < D; d++)
        pts_out (i, d) = center (d) + box_length * pts_in (i, d);
      wts_out (i) = box_factor * wts_in (i);
    }
}

/**
 * The maximum number of iterations for the Newton method.
 * If the Newton method does not converge within this number of iterations,
 * it is considered to have failed.
 */
const int NEWTON_ITER_TRESHOLD = 20;
/**
 * Maximum distance between linearly mapped point and Newton solution that is
 * considered acceptable (scaled w.r.t. element size).
 */
const double MAX_DIST_NEWTON = 5e-1;
/**
 * The tolerance for the Newton method.
 */
const double EPS_FIND_IP = 1e-12;

/**
 * Find the integration, i.e. the point in the reference element, that maps to
 * the given point in the physical element based on the given point in world
 * coordinates and the element transformation (constant or nonlinear).
 *
 * @tparam D The dimension of the vector.
 * @param ip The integration point to be found.
 * @param trafo The element transformation.
 * @param vec The vector describing the world coordinates.
 * @param lh The local heap.
 */
template <int D>
void FindIntegrationPoint (IntegrationPoint &ip,
                           const ElementTransformation &trafo,
                           const Vec<D> &vec, LocalHeap &lh)
{
  // HeapReset hr (lh);
  ip.Point ().Range (0, D) = 0.0;

  FlatVector<> diff (D, lh);
  FlatVector<> update (D, lh);

  // Linear guess
  auto mip_x0 = new (lh) MappedIntegrationPoint<D, D> (ip, trafo);
  diff = vec - mip_x0->GetPoint ();
  update = mip_x0->GetJacobianInverse () * diff;
  ip.Point ().Range (0, D) += update;
  const double h = sqrt (mip_x0->GetJacobiDet ());

  FlatVector<> ip_vec_lin (D, lh);
  ip_vec_lin = ip.Point ().Range (0, D);

  int its = 0;

  while (its == 0
         || (L2Norm (diff) > EPS_FIND_IP * h && its < NEWTON_ITER_TRESHOLD))
    {
      HeapReset hr (lh);
      auto mip_x0 = new (lh) MappedIntegrationPoint<D, D> (ip, trafo);
      diff = vec - mip_x0->GetPoint ();
      update = mip_x0->GetJacobianInverse () * diff;
      ip.Point ().Range (0, D) += update;
      its++;
    }

  FlatVector<> ip_vec (D, lh);
  ip_vec = ip.Point ().Range (0, D);

  if (its >= NEWTON_ITER_TRESHOLD || L2Norm (diff) > EPS_FIND_IP * h)
    throw Exception ("Newton did not converge after " + to_string (its)
                     + " iterations! (" + to_string (D) + "D)");
  else if (L2Norm (ip_vec - ip_vec_lin) > MAX_DIST_NEWTON * h)
    {
      cout << IM (6) << "Distance warning triggered, dist = "
           << L2Norm (ip_vec - ip_vec_lin) << " its = " << its << endl;
      cout << IM (6) << "taking a low order guess" << endl;
      ip.Point ().Range (0, D) = ip_vec_lin;
    }
  ip.SetFacetNr (-1, VOL);
}

/**
 * \brief Computes the points and weights for numerical integration over a box
 * part of an element or over the boundary of the box part.
 *
 * \tparam D The dimension of the (box) element.
 * \param order The order of the integration rule.
 * \param lh The local heap for temporary memory allocation.
 * \param element_vb The type of the box integration rule (VOL or BND).
 * \return A tuple containing the matrix of reference points and the vector of
 * reference weights.
 *
 * This function computes the points and weights for numerical integration over
 * a box element. The integration is performed using an integration rule of the
 * specified order. The points and weights are computed separately for the
 * volume (VOL) and boundary (BND) types of elements. The computed
 * points/weights are relative to the center of the box, with the box size
 * being provided.
 */
template <int D>
tuple<FlatMatrixFixWidth<D>, FlatVector<>>
GetBoxPointsAndWeights (int order, LocalHeap &lh, VorB element_vb = VOL,
                        BOXTYPE type = DEFAULT)
{
  const IntegrationRule &ir1d = SelectIntegrationRule (ET_SEGM, order);
  const int nip1D = ir1d.Size ();

  if ((type == BOXTYPE::BALL) && (D == 2))
    {
      const int ballnip = (order + 1) / 2;
      const int nip = nip1D * ballnip;
      FlatMatrixFixWidth<D> ref_pts_mat (nip, lh);
      FlatVector<> ref_wts_vec (nip, lh);
      for (int k = 0; k < ballnip; k++)
        {
          double etak = cos ((k + 1) * M_PI / (ballnip + 1));
          double ak
              = M_PI / (ballnip + 1) * sin ((k + 1) * M_PI / (ballnip + 1));
          for (int i : Range (nip1D))
            {
              int idx = nip1D * k + i;

              double segp = sqrt (1.0 - etak * etak);
              double p0 = ref_pts_mat (idx, 0) = etak;
              double p1 = ref_pts_mat (idx, 1)
                  = ir1d[i](0) * 2.0 * segp - segp;

              ref_wts_vec (idx) = pow (1.0 - (p0 * p0 + p1 * p1), 3) * ak
                                  * ir1d[i].Weight () * 2.0 * segp;
            }
        }
      return make_tuple (ref_pts_mat, ref_wts_vec);
    }

  switch (element_vb)
    {
    case VOL:
      {
        const int nip = pow (nip1D, D);
        FlatMatrixFixWidth<D> ref_pts_mat (nip, lh);
        FlatVector<> ref_wts_vec (nip, lh);
        for (int i : Range (nip))
          {
            int idx = i;
            ref_wts_vec (i) = 1.0;
            for (int d : Range (D))
              {
                ref_pts_mat (i, d) = ir1d[idx % nip1D](0) - 0.5;
                ref_wts_vec (i) *= ir1d[idx % nip1D].Weight ();
                idx = idx / nip1D;
              }
          }
        return make_tuple (ref_pts_mat, ref_wts_vec);
      }
      break;
    case BND:
      {
        const int nip_f = pow (nip1D, D - 1); // points per facet
        const int nip = 2 * D * nip_f;
        FlatMatrixFixWidth<D> ref_pts_mat (nip, lh);
        ref_pts_mat = 0.0;
        FlatVector<> ref_wts_vec (nip, lh);
        ref_wts_vec = 0.0;
        for (int bd : Range (D))
          {
            for (int i : Range (nip_f))
              {
                int idx = i;
                ref_wts_vec (2 * bd * nip_f + i) = 1.0;
                ref_pts_mat (2 * bd * nip_f + i, bd) = -0.5;
                ref_wts_vec (2 * bd * nip_f + i + nip_f) = 1.0;
                ref_pts_mat (2 * bd * nip_f + i + nip_f, bd) = 0.5;
                for (int d : Range (D))
                  {
                    if (d == bd)
                      continue;
                    ref_pts_mat (2 * bd * nip_f + i, d)
                        = ir1d[idx % nip1D](0) - 0.5;
                    ref_pts_mat (2 * bd * nip_f + i + nip_f, d)
                        = ir1d[idx % nip1D](0) - 0.5;
                    ref_wts_vec (2 * bd * nip_f + i)
                        *= ir1d[idx % nip1D].Weight ();
                    ref_wts_vec (2 * bd * nip_f + i + nip_f)
                        *= ir1d[idx % nip1D].Weight ();
                    idx = idx / nip1D;
                  }
              }
          }
        return make_tuple (ref_pts_mat, ref_wts_vec);
      }
      break;
    default:
      throw Exception ("GetBoxPointsAndWeights :: unhandled element_vb");
      break;
    }
  return make_tuple (FlatMatrixFixWidth<D> (), FlatVector<> ());
}

BoxIntegral ::BoxIntegral (shared_ptr<CoefficientFunction> _cf,
                           shared_ptr<BoxDifferentialSymbol> _dx)
    : Integral (_cf, *_dx), box_length (_dx->box_length),
      scale_with_elsize (_dx->scale_with_elsize), boxtype (_dx->boxtype)
{
  ;
}

BoxIntegral ::BoxIntegral (shared_ptr<CoefficientFunction> _cf,
                           DifferentialSymbol _dx, double _box_length,
                           bool _scale_with_elsize, BOXTYPE _boxtype)
    : Integral (_cf, _dx), box_length (_box_length),
      scale_with_elsize (_scale_with_elsize), boxtype (_boxtype)
{
  ;
}

shared_ptr<BilinearFormIntegrator>
BoxIntegral ::MakeBilinearFormIntegrator () const
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction &cf) {
    if (dynamic_cast<ProxyFunction *> (&cf))
      if (dynamic_cast<ProxyFunction &> (cf).IsOther ())
        has_other = true;
  });

  if (has_other)
    throw Exception ("no other terms in BoxIntegral..");
  if (dx.vb != VOL)
    throw Exception ("only VOL in BoxIntegral..");
  if (dx.skeleton)
    throw Exception ("no skeleton in BoxIntegral..");

  shared_ptr<BilinearFormIntegrator> bfi;
  bfi = make_shared<BoxBilinearFormIntegrator> (cf, dx.element_vb, box_length,
                                                scale_with_elsize, boxtype);

  if (dx.definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon);
          definedon_bitarray)
        bfi->SetDefinedOn (*definedon_bitarray);
    }
  bfi->SetDeformation (dx.deformation);
  bfi->SetBonusIntegrationOrder (dx.bonus_intorder);
  if (dx.definedonelements)
    bfi->SetDefinedOnElements (dx.definedonelements);
  return bfi;
}

shared_ptr<LinearFormIntegrator>
BoxIntegral ::MakeLinearFormIntegrator () const
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction &cf) {
    if (dynamic_cast<ProxyFunction *> (&cf))
      if (dynamic_cast<ProxyFunction &> (cf).IsOther ())
        has_other = true;
  });

  if (has_other)
    throw Exception ("no other terms in BoxIntegral..");
  if (dx.vb != VOL)
    throw Exception ("only VOL in BoxIntegral..");
  if (dx.skeleton)
    throw Exception ("no skeleton in BoxIntegral..");

  shared_ptr<LinearFormIntegrator> lfi;
  lfi = make_shared<BoxLinearFormIntegrator> (cf, dx.element_vb, box_length,
                                              scale_with_elsize, boxtype);
  if (dx.definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon);
          definedon_bitarray)
        lfi->SetDefinedOn (*definedon_bitarray);
    }
  lfi->SetDeformation (dx.deformation);
  lfi->SetBonusIntegrationOrder (dx.bonus_intorder);
  if (dx.definedonelements)
    lfi->SetDefinedOnElements (dx.definedonelements);

  return lfi;
}

template <typename TSCAL, int D>
TSCAL BoxIntegral ::T_BoxIntegrate (const ngcomp::MeshAccess &ma,
                                    FlatVector<TSCAL> element_wise)
{
  static Timer timer ("BoxIntegral::T_BoxIntegrate");
  RegionTimer reg (timer);
  LocalHeap glh (1000000000, "lh-T_BoxIntegrate");

  BitArray defon;
  if (dx.definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon))
        defon = *definedon_bitarray;
      if (auto definedon_string = get_if<string> (&*dx.definedon))
        {
          shared_ptr<MeshAccess> spma (const_cast<MeshAccess *> (&ma),
                                       NOOP_Deleter);
          Region reg (spma, dx.vb, *definedon_string);
          defon = reg.Mask ();
        }
    }

  int cfdim = cf->Dimension ();
  if (cfdim != 1)
    throw Exception (
        "only implemented for 1 dimensional coefficientfunctions");

  int order = 5 + dx.bonus_intorder;

  // auto [ref_pts_mat, ref_wts_vec] = GetBoxPointsAndWeights<D> (order, glh);
  tuple<FlatMatrixFixWidth<D>, FlatVector<>> paw
      = GetBoxPointsAndWeights<D> (order, glh, dx.element_vb, boxtype);
  auto ref_pts_mat = get<0> (paw);
  auto ref_wts_vec = get<1> (paw);
  const int nip = ref_pts_mat.Height ();

  try
    {
      TSCAL sum = 0.0;
      ma.IterateElements (VOL, glh, [&] (Ngs_Element el, LocalHeap &lh) {
        if (defon.Size () && !defon.Test (el.GetIndex ()))
          return;
        if (dx.definedonelements && !dx.definedonelements->Test (el.Nr ()))
          return;

        auto &trafo1 = ma.GetTrafo (el, lh);
        auto &trafo = trafo1.AddDeformation (this->dx.deformation.get (), lh);

        FlatMatrixFixWidth<D> pts_mat (nip, lh);
        FlatVector<> wts_vec (nip, lh);
        MapRefPoints<D> (trafo, ref_pts_mat, pts_mat, ref_wts_vec, wts_vec,
                         box_length, scale_with_elsize, lh, dx.element_vb);

        IntegrationRule &ir = *(new (lh) IntegrationRule (nip, lh));
        for (int i = 0; i < nip; i++)
          {
            FindIntegrationPoint<D> (ir[i], trafo, pts_mat.Row (i), lh);
            auto &mip
                = *(new (lh) MappedIntegrationPoint<D, D> ((ir)[i], trafo));
            ir[i].SetWeight (wts_vec (i) / mip.GetMeasure ());
            ir[i].SetNr (i);
          }

        SIMD_IntegrationRule &simd_ir
            = *(new (lh) SIMD_IntegrationRule (ir, lh));
        {
          SIMD_BaseMappedIntegrationRule &simd_mir = trafo (simd_ir, lh);
          FlatMatrix<SIMD<TSCAL>> val (simd_mir.Size (), 1, lh);
          cf->Evaluate (simd_mir, val);
          SIMD<TSCAL> lsum (0.0);
          for (size_t i = 0; i < simd_mir.Size (); i++)
            lsum += simd_mir[i].GetWeight () * val (i, 0);
          if (element_wise.Size ())
            element_wise (el.Nr ()) += HSum (lsum); // problem?
          AtomicAdd (sum, HSum (lsum));
        }
      });
      return ma.GetCommunicator ().AllReduce (sum, NG_MPI_SUM);
    }
  catch (ExceptionNOSIMD const &e)
    {
      cout << IM (6) << e.What () << "switching to non-SIMD evaluation"
           << endl;
    }

  TSCAL sum = 0.0;
  ma.IterateElements (VOL, glh, [&] (Ngs_Element el, LocalHeap &lh) {
    if (defon.Size () && !defon.Test (el.GetIndex ()))
      return;
    if (dx.definedonelements && !dx.definedonelements->Test (el.Nr ()))
      return;

    auto &trafo1 = ma.GetTrafo (el, lh);
    auto &trafo = trafo1.AddDeformation (this->dx.deformation.get (), lh);

    FlatMatrixFixWidth<D> pts_mat (nip, lh);
    FlatVector<> wts_vec (nip, lh);
    MapRefPoints<D> (trafo, ref_pts_mat, pts_mat, ref_wts_vec, wts_vec,
                     box_length, scale_with_elsize, lh, dx.element_vb);

    IntegrationRule &ir = *(new (lh) IntegrationRule (nip, lh));
    TSCAL lsum (0.0);
    for (int i = 0; i < nip; i++)
      {
        FindIntegrationPoint<D> (ir[i], trafo, pts_mat.Row (i), lh);
        ir[i].SetNr (i);
        auto &mip = *(new (lh) MappedIntegrationPoint<D, D> (ir[i], trafo));
        ir[i].SetWeight (wts_vec (i) / mip.GetMeasure ());
        auto &mip2 = *(new (lh) MappedIntegrationPoint<D, D> (ir[i], trafo));
        lsum += wts_vec (i) * cf->Evaluate (mip2);
      }
    if (element_wise.Size ())
      element_wise (el.Nr ()) += lsum;
    AtomicAdd (sum, lsum);
  });
  return ma.GetCommunicator ().AllReduce (sum, NG_MPI_SUM);
}

double BoxIntegral::Integrate (const ngcomp::MeshAccess &ma,
                               FlatVector<double> element_wise)
{
  double ret = 0.0;
  Switch<3> (ma.GetDimension () - 1, [&] (auto Dm1) {
    ret = T_BoxIntegrate<double, Dm1 + 1> (ma, element_wise);
  });
  return ret;
}

Complex BoxIntegral::Integrate (const ngcomp::MeshAccess &ma,
                                FlatVector<Complex> element_wise)
{
  Complex ret = Complex (0.0);
  Switch<3> (ma.GetDimension () - 1, [&] (auto Dm1) {
    ret = T_BoxIntegrate<Complex, Dm1 + 1> (ma, element_wise);
  });
  return ret;
}

template double
BoxIntegral ::T_BoxIntegrate<double, 1> (const ngcomp::MeshAccess &ma,
                                         FlatVector<double> element_wise);
template Complex
BoxIntegral ::T_BoxIntegrate<Complex, 1> (const ngcomp::MeshAccess &ma,
                                          FlatVector<Complex> element_wise);
template double
BoxIntegral ::T_BoxIntegrate<double, 2> (const ngcomp::MeshAccess &ma,
                                         FlatVector<double> element_wise);
template Complex
BoxIntegral ::T_BoxIntegrate<Complex, 2> (const ngcomp::MeshAccess &ma,
                                          FlatVector<Complex> element_wise);
template double
BoxIntegral ::T_BoxIntegrate<double, 3> (const ngcomp::MeshAccess &ma,
                                         FlatVector<double> element_wise);
template Complex
BoxIntegral ::T_BoxIntegrate<Complex, 3> (const ngcomp::MeshAccess &ma,
                                          FlatVector<Complex> element_wise);

BoxLinearFormIntegrator ::BoxLinearFormIntegrator (
    shared_ptr<CoefficientFunction> acf, VorB vb, double _box_length,
    bool _scale_with_elsize, BOXTYPE _boxtype)
    : SymbolicLinearFormIntegrator (acf, VOL, vb), box_length (_box_length),
      scale_with_elsize (_scale_with_elsize), boxtype (_boxtype)
{
  ;
}

void BoxLinearFormIntegrator ::CalcElementVector (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatVector<double> elvec, LocalHeap &lh) const
{
  Switch<3> (fel.Dim () - 1, [&] (auto Dm1) {
    T_CalcElementVector<Dm1 + 1> (fel, trafo, elvec, lh);
  });
}

void BoxLinearFormIntegrator ::CalcElementVector (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatVector<Complex> elvec, LocalHeap &lh) const
{
  Switch<3> (fel.Dim () - 1, [&] (auto Dm1) {
    T_CalcElementVector<Dm1 + 1> (fel, trafo, elvec, lh);
  });
}

template <int D, typename SCAL>
void BoxLinearFormIntegrator ::T_CalcElementVector (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatVector<SCAL> elvec, LocalHeap &lh) const
{
  static Timer timer ("BoxLFI - CalcElementVector");
  RegionTimer reg (timer);
  HeapReset hr (lh);

  int intorder = 2 * fel.Order () + bonus_intorder;

  tuple<FlatMatrixFixWidth<D>, FlatVector<>> paw
      = GetBoxPointsAndWeights<D> (intorder, lh, element_vb, boxtype);
  auto ref_pts_mat = get<0> (paw);
  auto ref_wts_vec = get<1> (paw);
  const int nip = ref_pts_mat.Height ();

  FlatMatrixFixWidth<D> pts_mat (nip, lh);
  FlatVector<> wts_vec (nip, lh);
  MapRefPoints<D> (trafo, ref_pts_mat, pts_mat, ref_wts_vec, wts_vec,
                   box_length, scale_with_elsize, lh, element_vb);

  IntegrationRule *ir = new (lh) IntegrationRule (nip, lh);
  for (int i = 0; i < nip; i++)
    {
      FindIntegrationPoint<D> ((*ir)[i], trafo, pts_mat.Row (i), lh);
      (*ir)[i].SetNr (i);
      auto &mip = *(new (lh) MappedIntegrationPoint<D, D> ((*ir)[i], trafo));
      (*ir)[i].SetWeight (wts_vec (i) / mip.GetMeasure ());
    }

  ProxyUserData ud (proxies.Size (), gridfunction_cfs.Size (), lh);
  const_cast<ElementTransformation &> (trafo).userdata = &ud;
  ud.fel = &fel;

  elvec = 0;
  if (simd_evaluate)
    {
      try
        {

          SIMD_IntegrationRule &simd_ir
              = *(new (lh) SIMD_IntegrationRule (*ir, lh));
          auto &simd_mir = trafo (simd_ir, lh);

          PrecomputeCacheCF (cache_cfs, simd_mir, lh);
          for (CoefficientFunction *cf : gridfunction_cfs)
            ud.AssignMemory (cf, simd_ir.GetNIP (), cf->Dimension (), lh);

          for (auto proxy : proxies)
            {
              FlatMatrix<SIMD<SCAL>> proxyvalues (proxy->Dimension (),
                                                  simd_ir.Size (), lh);
              for (size_t k = 0; k < proxy->Dimension (); k++)
                {
                  ud.testfunction = proxy;
                  ud.test_comp = k;

                  cf->Evaluate (simd_mir, proxyvalues.Rows (k, k + 1));
                  for (size_t i = 0; i < simd_mir.Size (); i++)
                    proxyvalues (k, i) *= simd_mir[i].GetWeight ();
                }

              proxy->Evaluator ()->AddTrans (fel, simd_mir, proxyvalues,
                                             elvec);
            }
        }
      catch (ExceptionNOSIMD const &e)
        {
          cout << IM (6) << e.What () << endl
               << "switching back to standard evaluation" << endl;
          simd_evaluate = false;
          T_CalcElementVector<D> (fel, trafo, elvec, lh);
        }
      return;
    }

  BaseMappedIntegrationRule &mir = trafo (*ir, lh);

  PrecomputeCacheCF (cache_cfs, mir, lh);
  for (CoefficientFunction *cf : gridfunction_cfs)
    ud.AssignMemory (cf, ir->GetNIP (), cf->Dimension (), lh);

  FlatVector<SCAL> elvec1 (elvec.Size (), lh);
  elvec1 = 0.0;
  FlatMatrix<SCAL> values (ir->Size (), 1, lh);

  for (auto proxy : proxies)
    {
      FlatMatrix<SCAL> proxyvalues (mir.Size (), proxy->Dimension (), lh);
      for (size_t k = 0; k < proxy->Dimension (); k++)
        {
          ud.testfunction = proxy;
          ud.test_comp = k;
          cf->Evaluate (mir, values);
          for (size_t i = 0; i < mir.Size (); i++)
            proxyvalues (i, k) = mir[i].GetWeight () * values (i, 0);
        }
      proxy->Evaluator ()->ApplyTrans (fel, mir, proxyvalues, elvec1, lh);
      elvec += elvec1;
    }

  return;
}

BoxBilinearFormIntegrator ::BoxBilinearFormIntegrator (
    shared_ptr<CoefficientFunction> acf, VorB avb, double _box_length,
    bool _scale_with_elsize, BOXTYPE _boxtype)
    : SymbolicBilinearFormIntegrator (acf, VOL, avb), box_length (_box_length),
      scale_with_elsize (_scale_with_elsize), boxtype (_boxtype)
{
  ;
}

void BoxBilinearFormIntegrator ::CalcElementMatrix (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatMatrix<double> elmat, LocalHeap &lh) const
{
  elmat = 0.0;
  Switch<3> (fel.Dim () - 1, [&] (auto Dm1) {
    T_CalcElementMatrixAdd<Dm1 + 1, double, double> (fel, trafo, elmat, lh);
  });
}

void BoxBilinearFormIntegrator ::CalcElementMatrixAdd (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatMatrix<double> elmat, bool &symmetric_so_far, LocalHeap &lh) const
{
  symmetric_so_far = false;
  Switch<3> (fel.Dim () - 1, [&] (auto Dm1) {
    T_CalcElementMatrixAdd<Dm1 + 1, double, double, double> (fel, trafo, elmat,
                                                             lh);
  });
}

void BoxBilinearFormIntegrator ::CalcElementMatrixAdd (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatMatrix<Complex> elmat, bool &symmetric_so_far, LocalHeap &lh) const
{
  symmetric_so_far = false;
  Switch<3> (fel.Dim () - 1, [&] (auto Dm1) {
    if (fel.ComplexShapes () || trafo.IsComplex ())
      T_CalcElementMatrixAdd<Dm1 + 1, Complex, Complex> (fel, trafo, elmat,
                                                         lh);
    else
      T_CalcElementMatrixAdd<Dm1 + 1, Complex, double> (fel, trafo, elmat, lh);
  });
}

template <int D, typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
void BoxBilinearFormIntegrator ::T_CalcElementMatrixAdd (
    const FiniteElement &fel, const ElementTransformation &trafo,
    FlatMatrix<SCAL_RES> elmat, LocalHeap &lh) const

{
  static Timer timer ("BoxBFI::CalcElementMatrixAdd");
  bool symmetric_so_far = false;

  auto save_userdata = trafo.PushUserData ();

  if (element_vb != VOL)
    {
      // T_CalcElementMatrixEBAdd<SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo,
      // elmat, lh); if (!IsSymmetric().IsTrue()) symmetric_so_far = false;
      // return;
      throw Exception ("only VOL in BoxIntegral..");
    }

  if (has_interpolate)
    {
      // T_CalcElementMatrixAddShapeWise<SCAL, SCAL_SHAPES, SCAL_RES> (fel,
      // trafo, elmat, lh); if (!IsSymmetric().IsTrue()) symmetric_so_far =
      // false; return;
      throw Exception ("no interpolate in BoxIntegral..");
    }

  bool is_mixedfe = typeid (fel) == typeid (const MixedFiniteElement &);
  const MixedFiniteElement *mixedfe
      = static_cast<const MixedFiniteElement *> (&fel);
  const FiniteElement &fel_trial = is_mixedfe ? mixedfe->FETrial () : fel;
  const FiniteElement &fel_test = is_mixedfe ? mixedfe->FETest () : fel;
  // size_t first_std_eval = 0;

  int intorder = fel_trial.Order () + fel_test.Order () + bonus_intorder;
  tuple<FlatMatrixFixWidth<D>, FlatVector<>> paw
      = GetBoxPointsAndWeights<D> (intorder, lh, element_vb, boxtype);
  auto ref_pts_mat = get<0> (paw);
  auto ref_wts_vec = get<1> (paw);
  const int nip = ref_pts_mat.Height ();

  FlatMatrixFixWidth<D> pts_mat (nip, lh);
  FlatVector<> wts_vec (nip, lh);
  MapRefPoints<D> (trafo, ref_pts_mat, pts_mat, ref_wts_vec, wts_vec,
                   box_length, scale_with_elsize, lh, element_vb);

  IntegrationRule *iir = new (lh) IntegrationRule (nip, lh);
  for (int i = 0; i < nip; i++)
    {
      FindIntegrationPoint<D> ((*iir)[i], trafo, pts_mat.Row (i), lh);
      auto &mip = *(new (lh) MappedIntegrationPoint<D, D> ((*iir)[i], trafo));
      (*iir)[i].SetWeight (wts_vec (i) / mip.GetMeasure ());
      (*iir)[i].SetNr (i);
    }

  if (simd_evaluate)
    try
      {
        SIMD_IntegrationRule &ir = *(new (lh) SIMD_IntegrationRule (*iir, lh));
        SIMD_BaseMappedIntegrationRule &mir = trafo (ir, lh);

        ProxyUserData ud;
        const_cast<ElementTransformation &> (trafo).userdata = &ud;
        PrecomputeCacheCF (cache_cfs, mir, lh);

        // bool symmetric_so_far = true;
        int k1 = 0;
        int k1nr = 0;

        for (auto proxy1 : trial_proxies)
          {
            int l1 = 0;
            int l1nr = 0;
            for (auto proxy2 : test_proxies)
              {
                size_t dim_proxy1 = proxy1->Dimension ();
                size_t dim_proxy2 = proxy2->Dimension ();

                size_t tt_pair = l1nr * trial_proxies.Size () + k1nr;
                // first_std_eval = k1nr*test_proxies.Size()+l1nr;  // in case
                // of SIMDException
                bool is_nonzero = nonzeros_proxies (tt_pair);
                bool is_diagonal = diagonal_proxies (tt_pair);

                if (is_nonzero)
                  {
                    HeapReset hr (lh);
                    bool samediffop = same_diffops (tt_pair) && !is_mixedfe;
                    // td.Start();

                    FlatMatrix<SIMD<SCAL>> proxyvalues (
                        dim_proxy1 * dim_proxy2, ir.Size (), lh);
                    FlatMatrix<SIMD<SCAL>> diagproxyvalues (dim_proxy1,
                                                            ir.Size (), lh);
                    FlatMatrix<SIMD<SCAL>> val (1, ir.Size (), lh);

                    IntRange r1 = proxy1->Evaluator ()->UsedDofs (fel_trial);
                    IntRange r2 = proxy2->Evaluator ()->UsedDofs (fel_test);
                    SliceMatrix<SCAL_RES> part_elmat
                        = elmat.Rows (r2).Cols (r1);

                    FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1 (
                        elmat.Width () * dim_proxy1, ir.Size (), lh);
                    FlatMatrix<SIMD<SCAL>> bdbmat1 (
                        elmat.Width () * dim_proxy2, ir.Size (), lh);
                    FlatMatrix<SIMD<SCAL_SHAPES>> bbmat2
                        = samediffop ? bbmat1
                                     : FlatMatrix<SIMD<SCAL_SHAPES>> (
                                           elmat.Height () * dim_proxy2,
                                           ir.Size (), lh);

                    FlatMatrix<SIMD<SCAL>> hbdbmat1 (elmat.Width (),
                                                     dim_proxy2 * ir.Size (),
                                                     bdbmat1.Data ());
                    FlatMatrix<SIMD<SCAL_SHAPES>> hbbmat2 (
                        elmat.Height (), dim_proxy2 * ir.Size (),
                        bbmat2.Data ());

                    if (ddcf_dtest_dtrial (l1nr, k1nr))
                      {
                        ddcf_dtest_dtrial (l1nr, k1nr)
                            ->Evaluate (mir, proxyvalues);

                        if (is_diagonal)
                          for (auto k : Range (dim_proxy1))
                            diagproxyvalues.Row (k)
                                = proxyvalues.Row (k * (dim_proxy1 + 1));
                      }
                    else
                      {
                        //                          RegionTimer regdmat(tdmat);
                        if (!is_diagonal)
                          {
                            for (size_t k = 0, kk = 0; k < dim_proxy1; k++)
                              for (size_t l = 0; l < dim_proxy2; l++, kk++)
                                {
                                  if (nonzeros (l1 + l, k1 + k))
                                    {
                                      ud.trialfunction = proxy1;
                                      ud.trial_comp = k;
                                      ud.testfunction = proxy2;
                                      ud.test_comp = l;

                                      cf->Evaluate (
                                          mir, proxyvalues.Rows (kk, kk + 1));
                                    }
                                  else
                                    {
                                      ;
                                    }
                                  // proxyvalues.Row(kk) = 0.0;
                                }
                          }
                        else
                          {
                            for (size_t k = 0; k < dim_proxy1; k++)
                              {
                                ud.trialfunction = proxy1;
                                ud.trial_comp = k;
                                ud.testfunction = proxy2;
                                ud.test_comp = k;

                                cf->Evaluate (mir,
                                              diagproxyvalues.Rows (k, k + 1));
                              }
                          }
                        // td.Stop();
                      }

                    // NgProfiler::StartThreadTimer (timer_SymbBFIscale,
                    // TaskManager::GetThreadId());
                    FlatVector<SIMD<double>> weights (ir.Size (), lh);
                    if (!is_diagonal)
                      for (size_t i = 0; i < ir.Size (); i++)
                        // proxyvalues.Col(i) *= mir[i].GetWeight();
                        weights (i) = mir[i].GetWeight ();
                    else
                      for (size_t i = 0; i < ir.Size (); i++)
                        diagproxyvalues.Col (i) *= mir[i].GetWeight ();

                    {
                      // RegionTimer regbmat(timer_SymbBFIbmat);
                      proxy1->Evaluator ()->CalcMatrix (fel_trial, mir,
                                                        bbmat1);
                      if (!samediffop)
                        proxy2->Evaluator ()->CalcMatrix (fel_test, mir,
                                                          bbmat2);
                    }

                    if (is_diagonal)
                      {
                        // size_t sr1 = r1.Size();
                        for (size_t j = 0; j < dim_proxy1; j++)
                          {
                            auto hbbmat1
                                = bbmat1.RowSlice (j, dim_proxy1).Rows (r1);
                            auto hbdbmat1
                                = bdbmat1.RowSlice (j, dim_proxy1).Rows (r1);

                            for (size_t k = 0; k < bdbmat1.Width (); k++)
                              hbdbmat1.Col (k).Range (0, r1.Size ())
                                  = diagproxyvalues (j, k) * hbbmat1.Col (k);
                          }
                      }
                    else
                      {
                        hbdbmat1.Rows (r1) = 0.0;
                        for (size_t j = 0; j < dim_proxy2; j++)
                          for (size_t k = 0; k < dim_proxy1; k++)
                            if (nonzeros (l1 + j, k1 + k))
                              {
                                auto proxyvalues_jk
                                    = symbolic_integrator_uses_diff
                                          ? proxyvalues.Row (j * dim_proxy1
                                                             + k)
                                          : proxyvalues.Row (k * dim_proxy2
                                                             + j);
                                auto bbmat1_k = bbmat1.RowSlice (k, dim_proxy1)
                                                    .Rows (r1);
                                auto bdbmat1_j
                                    = bdbmat1.RowSlice (j, dim_proxy2)
                                          .Rows (r1);

                                for (size_t i = 0; i < ir.Size (); i++)
                                  bdbmat1_j.Col (i).Range (0, r1.Size ())
                                      += proxyvalues_jk (i) * weights (i)
                                         * bbmat1_k.Col (i);
                              }
                      }

                    symmetric_so_far &= samediffop && is_diagonal;

                    {
                      // static Timer t("AddABt", NoTracing);
                      // RegionTracer reg(TaskManager::GetThreadId(), t);

                      if (symmetric_so_far)
                        {
                          AddABtSym (hbbmat2.Rows (r2), hbdbmat1.Rows (r1),
                                     part_elmat);
                        }
                      else
                        {
                          AddABt (hbbmat2.Rows (r2), hbdbmat1.Rows (r1),
                                  part_elmat);
                        }
                    }
                    // if (symmetric_so_far)
                    //{
                    // ExtendSymmetric (part_elmat);
                    //}
                  }

                l1 += proxy2->Dimension ();
                l1nr++;
              }
            k1 += proxy1->Dimension ();
            k1nr++;
          }

        return;
      }
    catch (ExceptionNOSIMD const &e)
      {
        cout << IM (6) << e.What () << endl
             << "switching to scalar evaluation" << endl;
        simd_evaluate = false;
        T_CalcElementMatrixAdd<D, SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo,
                                                                elmat, lh);
        return;
      }

  BaseMappedIntegrationRule &mir = trafo (*iir, lh);

  ProxyUserData ud;
  const_cast<ElementTransformation &> (trafo).userdata = &ud;
  PrecomputeCacheCF (cache_cfs, mir, lh);

  // tstart.Stop();
  // bool symmetric_so_far = true;
  int k1 = 0;
  int k1nr = 0;
  for (auto proxy1 : trial_proxies)
    {
      int l1 = 0;
      int l1nr = 0;
      for (auto proxy2 : test_proxies)
        {
          bool is_diagonal = proxy1->Dimension () == proxy2->Dimension ();
          bool is_nonzero = false;

          for (size_t k = 0; k < proxy1->Dimension (); k++)
            for (size_t l = 0; l < proxy2->Dimension (); l++)
              if (nonzeros (l1 + l, k1 + k))
                {
                  if (k != l)
                    is_diagonal = false;
                  is_nonzero = true;
                }

          if (is_nonzero) //   && k1nr*test_proxies.Size()+l1nr >=
                          //   first_std_eval)
            {
              HeapReset hr (lh);
              bool samediffop
                  = (*(proxy1->Evaluator ()) == *(proxy2->Evaluator ()))
                    && !is_mixedfe;
              // td.Start();
              FlatTensor<3, SCAL> proxyvalues (
                  lh, mir.Size (), proxy1->Dimension (), proxy2->Dimension ());
              FlatVector<SCAL> diagproxyvalues (
                  mir.Size () * proxy1->Dimension (), lh);
              FlatMatrix<SCAL> val (mir.Size (), 1, lh);

              IntRange r1 = proxy1->Evaluator ()->UsedDofs (fel_trial);
              IntRange r2 = proxy2->Evaluator ()->UsedDofs (fel_test);
              SliceMatrix<SCAL_RES> part_elmat = elmat.Rows (r2).Cols (r1);
              FlatMatrix<SCAL_SHAPES, ColMajor> bmat1 (proxy1->Dimension (),
                                                       elmat.Width (), lh);
              FlatMatrix<SCAL_SHAPES, ColMajor> bmat2 (proxy2->Dimension (),
                                                       elmat.Height (), lh);

              if (ddcf_dtest_dtrial (l1nr, k1nr))
                {
                  //                    cout << "use ddcf_dtest_dtrial (NO
                  //                    SIMD)" << endl;
                  // TODO: optimize for element-wise constant case?
                  FlatMatrix<SCAL> mproxyvalues (
                      mir.Size (), proxy1->Dimension () * proxy2->Dimension (),
                      proxyvalues.Data ());
                  ddcf_dtest_dtrial (l1nr, k1nr)->Evaluate (mir, mproxyvalues);
                  if (is_diagonal)
                    for (auto k : Range (proxy1->Dimension ()))
                      diagproxyvalues.Slice (k, proxy1->Dimension ())
                          = proxyvalues (STAR, k, k);
                }
              else
                {
                  if (!is_diagonal)
                    for (size_t k = 0; k < proxy1->Dimension (); k++)
                      for (size_t l = 0; l < proxy2->Dimension (); l++)
                        {
                          if (nonzeros (l1 + l, k1 + k))
                            {
                              ud.trialfunction = proxy1;
                              ud.trial_comp = k;
                              ud.testfunction = proxy2;
                              ud.test_comp = l;

                              cf->Evaluate (mir, val);
                              proxyvalues (STAR, k, l) = val.Col (0);
                            }
                          else
                            proxyvalues (STAR, k, l) = 0.0;
                        }
                  else
                    for (size_t k = 0; k < proxy1->Dimension (); k++)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = k;

                        if (!elementwise_constant)
                          {
                            cf->Evaluate (mir, val);
                            diagproxyvalues.Slice (k, proxy1->Dimension ())
                                = val.Col (0);
                          }
                        else
                          {
                            cf->Evaluate (mir[0], val.Row (0));
                            diagproxyvalues.Slice (k, proxy1->Dimension ())
                                = val (0, 0);
                          }
                      }
                }
              // td.Stop();

              if (!mir.IsComplex ())
                {
                  if (!is_diagonal)
                    for (size_t i = 0; i < mir.Size (); i++)
                      proxyvalues (i, STAR, STAR) *= mir[i].GetWeight ();
                  else
                    for (size_t i = 0; i < mir.Size (); i++)
                      diagproxyvalues.Range (proxy1->Dimension ()
                                             * IntRange (i, i + 1))
                          *= mir[i].GetWeight ();
                }
              else
                { // pml
                  if (!is_diagonal)
                    for (size_t i = 0; i < mir.Size (); i++)
                      proxyvalues (i, STAR, STAR) *= mir[i].GetWeight ();
                  else
                    for (size_t i = 0; i < mir.Size (); i++)
                      diagproxyvalues.Range (proxy1->Dimension ()
                                             * IntRange (i, i + 1))
                          *= static_cast<
                                 const ScalMappedIntegrationPoint<SCAL> &> (
                                 mir[i])
                                 .GetJacobiDet ()
                             * (*iir)[i].Weight ();
                }

              constexpr size_t BS = 16;
              for (size_t i = 0; i < mir.Size (); i += BS)
                {
                  HeapReset hr (lh);
                  int bs = min2 (size_t (BS), mir.Size () - i);

                  FlatMatrix<SCAL_SHAPES> bbmat1 (
                      elmat.Width (), bs * proxy1->Dimension (), lh);
                  FlatMatrix<SCAL> bdbmat1 (elmat.Width (),
                                            bs * proxy2->Dimension (), lh);
                  FlatMatrix<SCAL_SHAPES> bbmat2
                      = samediffop ? bbmat1
                                   : FlatMatrix<SCAL_SHAPES> (
                                         elmat.Height (),
                                         bs * proxy2->Dimension (), lh);

                  // tb.Start();
                  BaseMappedIntegrationRule &bmir = mir.Range (i, i + bs, lh);

                  proxy1->Evaluator ()->CalcMatrix (fel_trial, bmir,
                                                    Trans (bbmat1), lh);

                  if (!samediffop)
                    proxy2->Evaluator ()->CalcMatrix (fel_test, bmir,
                                                      Trans (bbmat2), lh);
                  // tb.Stop();

                  // tdb.Start();
                  if (is_diagonal)
                    {
                      FlatVector<SCAL> diagd (bs * proxy1->Dimension (), lh);
                      diagd = diagproxyvalues.Range (
                          i * proxy1->Dimension (),
                          (i + bs) * proxy1->Dimension ());
                      for (size_t i = 0; i < diagd.Size (); i++)
                        bdbmat1.Col (i) = diagd (i) * bbmat1.Col (i);
                      // MultMatDiagMat(bbmat1, diagd, bdbmat1);
                      // tdb.AddFlops (bbmat1.Height()*bbmat1.Width());
                    }
                  else
                    {
                      for (int j = 0; j < bs; j++)
                        {
                          int ii = i + j;
                          IntRange r1
                              = proxy1->Dimension () * IntRange (j, j + 1);
                          IntRange r2
                              = proxy2->Dimension () * IntRange (j, j + 1);
                          // bdbmat1.Cols(r2) = bbmat1.Cols(r1) *
                          // proxyvalues(ii,STAR,STAR);
                          MultMatMat (bbmat1.Cols (r1),
                                      proxyvalues (ii, STAR, STAR),
                                      bdbmat1.Cols (r2));
                        }
                      // tdb.AddFlops
                      // (proxy1->Dimension()*proxy2->Dimension()*bs*bbmat1.Height());
                    }
                  // tdb.Stop();
                  // tlapack.Start();
                  // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) *
                  // Trans(bdbmat1.Rows(r1)); AddABt (bbmat2.Rows(r2),
                  // bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                  symmetric_so_far &= samediffop && is_diagonal;
                  if (symmetric_so_far)
                    AddABtSym (bbmat2.Rows (r2), bdbmat1.Rows (r1),
                               part_elmat);
                  else
                    AddABt (bbmat2.Rows (r2), bdbmat1.Rows (r1), part_elmat);
                  // tlapack.Stop();
                  // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                }

              if (symmetric_so_far)
                for (size_t i = 0; i < part_elmat.Height (); i++)
                  for (size_t j = i + 1; j < part_elmat.Width (); j++)
                    part_elmat (i, j) = part_elmat (j, i);
            }

          l1 += proxy2->Dimension ();
          l1nr++;
        }
      k1 += proxy1->Dimension ();
      k1nr++;
    }
}

void BoxBilinearFormIntegrator ::ApplyElementMatrix (
    const FiniteElement &, const ElementTransformation &,
    const FlatVector<double>, FlatVector<double>, void *, LocalHeap &) const
{
  throw Exception (
      "BoxBilinearFormIntegrator::ApplyElementMatrix not implemented");
}

void BoxBilinearFormIntegrator ::CalcLinearizedElementMatrix (
    const FiniteElement &, const ElementTransformation &, FlatVector<double>,
    FlatMatrix<double>, LocalHeap &) const
{
  throw Exception ("BoxBilinearFormIntegrator::CalcLinearizedElementMatrix "
                   "not implemented");
}

////////////////////////// python interface ///////////////////////////

#ifdef NGS_PYTHON

void ExportBoxIntegral (py::module m)
{

  py::enum_<BOXTYPE> (m, "BOXTYPE", docu_string (R"raw_string(
  Shape of subdomain for BoxIntegral, currently supported are:
  BOX: gives square or cube, in 2D or 3D
  BALL: gives circle, currently only in 2D
  )raw_string"))
      .value ("DEFAULT", DEFAULT)
      .value ("BOX", BOX)
      .value ("BALL", BALL);

  const auto py_boxint
      = py::class_<BoxIntegral, shared_ptr<BoxIntegral>, Integral> (
          m, "BoxIntegral", docu_string (R"raw_string(
        BoxIntegral allows to formulate linear, bilinear forms and integrals on
        box parts of the mesh")raw_string"));

  py::class_<BoxDifferentialSymbol, DifferentialSymbol> (
      m, "BoxDifferentialSymbol", docu_string (R"raw_string(
dBox that allows to formulate linear, bilinear forms and integrals on
(bounding) boxes

Example use case:

  dbox = BoxDifferentialSymbol()
  a = BilinearForm(...)
  a += u * v * dbox(element_boundary=...)

)raw_string"))
      .def (py::init<> (), docu_string (R"raw_string(
Constructor of BoxDifferentialSymbol.

  Argument: none
)raw_string"))
      .def (
          "__call__",
          [] (BoxDifferentialSymbol &self,
              optional<variant<Region, string>> definedon,
              bool element_boundary, VorB element_vb,
              shared_ptr<GridFunction> deformation,
              shared_ptr<BitArray> definedonelements, int bonus_intorder,
              double box_length, bool scale_with_elsize, BOXTYPE boxtype) {
            if (element_boundary)
              element_vb = BND;
            BOXTYPE newboxtype
                = boxtype >= BOXTYPE::DEFAULT ? boxtype : self.boxtype;
            auto dx = BoxDifferentialSymbol (box_length, scale_with_elsize,
                                             newboxtype);
            dx.element_vb = element_vb;
            dx.bonus_intorder = bonus_intorder;
            if (definedon)
              {
                if (auto definedon_region = get_if<Region> (&*definedon);
                    definedon_region)
                  {
                    dx.definedon = definedon_region->Mask ();
                    dx.vb = VorB (*definedon_region);
                  }
                if (auto definedon_string = get_if<string> (&*definedon);
                    definedon_string)
                  dx.definedon = *definedon_string;
              }
            dx.deformation = deformation;
            dx.definedonelements = definedonelements;
            return dx;
          },
          py::arg ("definedon") = nullptr,
          py::arg ("element_boundary") = false, py::arg ("element_vb") = VOL,
          py::arg ("deformation") = nullptr,
          py::arg ("definedonelements") = nullptr,
          py::arg ("bonus_intorder") = 0, py::arg ("box_length") = 0.5,
          py::arg ("scale_with_elsize") = false,
          py::arg ("boxtype") = BOXTYPE::DEFAULT, docu_string (R"raw_string(
The call of a BoxDifferentialSymbol allows to specify what is needed to specify the
integration domain. It returns a new BoxDifferentialSymbol.

Parameters:

definedon (Region or Array) : specifies on which part of the mesh (in terms of regions)
  the current form shall be defined.
element_boundary (bool) : Does the integral take place on the boundary of an element-
  boundary?
element_vb (VOL/BND) : Where does the integral take place from point of view
  of an element (BBND/BBBND are not implemented).
deformation (GridFunction) : which mesh deformation shall be applied (default : None)
definedonelements (BitArray) : Set of elements or facets where the integral shall be
  defined.
bonus_intorder (int) : additional integration order for the integration rule (default: 0)
box_length (double) : length of the box (default: 0.5)
scale_with_elsize (bool) : if true, the box length is scaled with the size of the
  element (default: false)
boxtype (BOXTYPE) : shape of the box (default: BOX)
)raw_string"))
      .def_property (
          "element_vb",
          [] (BoxDifferentialSymbol &self) { return self.element_vb; },
          [] (BoxDifferentialSymbol &self, VorB element_vb) {
            self.element_vb = element_vb;
            return self.element_vb;
          },
          "box volume or box boundary integral on each (volume) element?");
}
#endif // NGS_PYTHON
