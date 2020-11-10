#include <comp.hpp> // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <python_comp.hpp>
#include <fem.hpp>
#include <multigrid.hpp>

#include "trefftzwavefe.hpp"
#include "qtrefftzwavefe.hpp"
#include "trefftzfespace.hpp"
#include "diffopmapped.hpp"

namespace ngcomp
{

  TrefftzFESpace ::TrefftzFESpace (shared_ptr<MeshAccess> ama,
                                   const Flags &flags)
      : FESpace (ama, flags)
  {
    type = "trefftzfespace";

    // cout << "======== Constructor of TrefftzFESpace =========" << endl;
    // cout << "Flags:" << endl << flags;

    fullD = ma->GetDimension ();
    D = fullD - 1;

    this->dgjumps = true;
    order = int (flags.GetNumFlag ("order", 3));
    c = flags.GetNumFlag ("wavespeed", 1);
    basistype = flags.GetNumFlag ("basistype", 0);
    useshift = flags.GetNumFlag ("useshift", 1);
    usescale = flags.GetNumFlag ("usescale", 1);
    useqt = flags.GetNumFlag ("useqt", 0);

    local_ndof = (BinCoeff (fullD - 1 + order, order)
                  + BinCoeff (fullD - 1 + order - 1, order - 1));
    nel = ma->GetNE ();
    ndof = local_ndof * nel;

    switch (fullD)
      {
      case 1:
        {
          // evaluator[VOL] =
          // make_shared<T_DifferentialOperator<DiffOpMapped<1>>>();
          // flux_evaluator[VOL] =
          // make_shared<T_DifferentialOperator<DiffOpMappedGradient<1>>>();
          // TrefftzWaveBasis<1>::getInstance().CreateTB(order, basistype);
          // break;
        }
      case 2:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<2>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<2>>> ();
          TrefftzWaveBasis<1>::getInstance ().CreateTB (order, basistype);
          break;
        }
      case 3:
        {
          evaluator[VOL]
              = make_shared<T_DifferentialOperator<DiffOpMapped<3>>> ();
          flux_evaluator[VOL] = make_shared<
              T_DifferentialOperator<DiffOpMappedGradient<3>>> ();
          TrefftzWaveBasis<2>::getInstance ().CreateTB (order, basistype);
          break;
        }
      }

    switch (fullD)
      {
      case 2:
        additional_evaluators.Set (
            "hesse",
            make_shared<T_DifferentialOperator<DiffOpMappedHesse<2>>> ());
        break;
      case 3:
        additional_evaluators.Set (
            "hesse",
            make_shared<T_DifferentialOperator<DiffOpMappedHesse<3>>> ());
        break;
      default:;
      }
  }

  void TrefftzFESpace ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize (0);
    if (!DefinedOn (ei) || ei.VB () != VOL)
      return;
    // int n_vert = ma->GetNV();		int n_edge = ma->GetNEdges();		int n_cell =
    // ma->GetNE(); Ngs_Element ngel = ma->GetElement (ei);
    for (int j = ei.Nr () * local_ndof; j < local_ndof * (ei.Nr () + 1); j++)
      {
        dnums.Append (j);
      }
    // cout << "GetDofNrs: ei.Nr() = " << ei.Nr() << " local_ndof:" <<
    // local_ndof << " ndof: " << ndof << " dnums: \n" << dnums << endl <<
    // "================================================" << endl ;
  }

  FiniteElement &TrefftzFESpace ::GetFE (ElementId ei, Allocator &alloc) const
  {

    // auto vertices_index = ma->GetElVertices(ei);
    // cout << "element vectice coord: \n"  <<
    // ma->GetPoint<3>(vertices_index[0]) << endl<<
    // ma->GetPoint<3>(vertices_index[1])
    // <<endl<<ma->GetPoint<3>(vertices_index[2])<<endl<<ma->GetPoint<3>(vertices_index[3])<<endl;

    Ngs_Element ngel = ma->GetElement (ei);
    ELEMENT_TYPE eltype = ngel.GetType ();

    if (ei.IsVolume ())
      {
        switch (ma->GetElType (ei))
          {
          case ET_SEGM:
            {
              // return *(new (alloc)
              // TrefftzWaveFE<1>(order,c,ElCenter<1>(ei),Adiam<1>(ei),ET_SEGM));
              // break;
            }
          case ET_QUAD:
          case ET_TRIG:
            {
              LocalHeap lh (1000 * 1000);
              const ELEMENT_TYPE eltyp = ET_TRIG;
              const int D = 2;
              IntegrationRule ir (eltyp, 0);
              MappedIntegrationPoint<D, D> mip (
                  ir[0], ma->GetTrafo (ElementId (0), lh));
              mip.Point () = ElCenter<1> (ei).Range (0, 1);
              if (useqt)
                {
                  // Matrix<> gamma(this->order,this->order);
                  // for(int nx=0;nx<this->order-1;nx++)
                  //{
                  // int ny = 0;
                  // gamma(nx,ny) =
                  // wavespeedmatrix(nx,ny)->Evaluate(mip)/(factorial(nx)*factorial(ny));
                  //}
                  return *(new (alloc) QTrefftzWaveFE<1> (
                      this->gamma[ei.Nr ()], order, ElCenter<1> (ei), 1.0,
                      ma->GetElType (ei)));
                }
              else
                return *(new (alloc) TrefftzWaveFE<1> (
                    order,
                    wavespeedcf != NULL ? wavespeedcf->Evaluate (mip) : c,
                    ElCenter<1> (ei),
                    Adiam<1> (ei, wavespeedcf != NULL
                                      ? wavespeedcf->Evaluate (mip)
                                      : c),
                    ma->GetElType (ei)));
              break;
            }
          case ET_HEX:
          case ET_PRISM:
          case ET_PYRAMID:
          case ET_TET:
            {
              LocalHeap lh (1000 * 1000);
              const ELEMENT_TYPE eltyp = ET_TRIG;
              const int D = 3;
              IntegrationRule ir (eltyp, 0);
              MappedIntegrationPoint<D, D> mip (
                  ir[0], ma->GetTrafo (ElementId (0), lh));
              mip.Point () = ElCenter<2> (ei).Range (0, 2);

              if (useqt)
                {
                  static Timer timereval ("evalc", 2);
                  timereval.Start ();
                  // Matrix<> gamma(this->order-1);
                  // for(int nx=0;nx<this->order-1;nx++)
                  //{
                  // for(int ny=0;ny<this->order-1;ny++)
                  //{
                  // gamma(nx,ny) =
                  // wavespeedmatrix(nx,ny)->Evaluate(mip)/(factorial(nx)*factorial(ny));
                  //}
                  //}
                  timereval.Stop ();
                  return *(new (alloc) QTrefftzWaveFE<2> (
                      this->gamma[ei.Nr ()], order, ElCenter<2> (ei), 1.0,
                      ma->GetElType (ei)));
                }
              else
                return *(new (alloc) TrefftzWaveFE<2> (
                    order,
                    wavespeedcf != NULL ? wavespeedcf->Evaluate (mip) : c,
                    ElCenter<2> (ei),
                    Adiam<2> (ei, wavespeedcf != NULL
                                      ? wavespeedcf->Evaluate (mip)
                                      : c),
                    ma->GetElType (ei)));
            }
            break;
          }
      }
    else
      {
        try
          {
            return SwitchET<ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD> (
                eltype, [&alloc] (auto et) -> FiniteElement & {
                  return *new (alloc) DummyFE<et.ElementType ()>;
                });
          }
        catch (Exception e)
          {
            throw Exception ("illegal element type in Trefftz::GetSurfaceFE");
          }
      }
  }

  template <int D> double TrefftzFESpace ::Adiam (ElementId ei, double c) const
  {
    double anisotropicdiam = 0.0;
    auto vertices_index = ma->GetElVertices (ei);
    for (auto vertex1 : vertices_index)
      {
        for (auto vertex2 : vertices_index)
          {
            Vec<D + 1> v1 = ma->GetPoint<D + 1> (vertex1);
            Vec<D + 1> v2 = ma->GetPoint<D + 1> (vertex2);
            // cout << "v1: " << v1 << " v1 part: " << v1(1,D-1) << "norm " <<
            // L2Norm(v1) << endl ;
            anisotropicdiam
                = max (anisotropicdiam,
                       sqrt (L2Norm2 (v1.Range (0, D) - v2.Range (0, D))
                             + pow (c * (v1 (D) - v2 (D)), 2)));
          }
      }
    return anisotropicdiam * usescale + (usescale == 0);
  }

  template <int D>
  double TrefftzFESpace ::Adiam (ElementId ei,
                                 shared_ptr<CoefficientFunction> c) const
  {
    LocalHeap lh (1000 * 1000);
    double anisotropicdiam = 0.0;
    auto vertices_index = ma->GetElVertices (ei);

    for (auto vertex1 : vertices_index)
      {
        for (auto vertex2 : vertices_index)
          {
            Vec<D + 1> v1 = ma->GetPoint<D + 1> (vertex1);
            Vec<D + 1> v2 = ma->GetPoint<D + 1> (vertex2);
            // cout << "v1: " << v1 << " v1 part: " << v1.Range(0,D) << " el
            // type: " << ma->GetElType(ei) << " norm " << L2Norm(v1) << endl ;
            IntegrationRule ir (ma->GetElType (ei), 0);
            ElementTransformation &trafo = ma->GetTrafo (ei, lh);
            MappedIntegrationPoint<D + 1, D + 1> mip (ir[0], trafo);
            mip.Point () = v1;
            double c1 = wavespeedcf->Evaluate (mip);
            mip.Point () = v2;
            double c2 = wavespeedcf->Evaluate (mip);

            anisotropicdiam
                = max (anisotropicdiam,
                       sqrt (L2Norm2 (v1.Range (0, D) - v2.Range (0, D))
                             + pow (c1 * v1 (D) - c2 * v2 (D), 2)));
          }
      }
    return anisotropicdiam * usescale + (usescale == 0);
  }

  template <int D> Vec<D + 1> TrefftzFESpace ::ElCenter (ElementId ei) const
  {
    Vec<D + 1> center = 0;
    auto vertices_index = ma->GetElVertices (ei);
    for (auto vertex : vertices_index)
      center += ma->GetPoint<D + 1> (vertex);
    center *= (1.0 / vertices_index.Size ()) * useshift;
    return center;
  }

  DocInfo TrefftzFESpace ::GetDocu ()
  {
    auto docu = FESpace::GetDocu ();
    docu.Arg ("useshift")
        = "bool = True\n"
          "  use shift of basis functins to element center and scale them";
    docu.Arg ("basistype")
        = "bool = True\n"
          "  use shift of basis functins to element center and scale them";
    docu.Arg ("wavespeed")
        = "bool = True\n"
          "  use shift of basis functins to element center and scale them";
    return docu;
  }

  /*
     register fe-spaces
     Object of type TrefftzFESpace can be defined in the pde-file via
     "define fespace v -type=trefftzfespace"
     */
  static RegisterFESpace<TrefftzFESpace> initi_trefftz ("trefftzfespace");

}

#ifdef NGS_PYTHON
void ExportTrefftzFESpace (py::module m)
{
  using namespace ngcomp;
  // using namespace ngfem;
  //[>
  // We just export the class here and use the FESpace constructor to create
  // our space. This has the advantage, that we do not need to specify all the
  // flags to parse (like dirichlet, definedon,...), but we can still append
  // new functions only for that space.
  //*/
  // py::class_<TrefftzFESpace, shared_ptr<TrefftzFESpace>, FESpace>
  //(m, "TrefftzFESpace", "FESpace with first order and second order trigs on
  //2d mesh") .def("GetNDof", &TrefftzFESpace::GetNDof)
  //;
  // m.def("GetNDof", [](shared_ptr<FESpace> fes) {
  // cout << typeid(*fes).name() << endl;
  ////fes->GetNDof();
  //});

  ExportFESpace<TrefftzFESpace> (m, "trefftzfespace")
      .def ("GetDocu", &TrefftzFESpace::GetDocu)
      .def ("GetNDof", &TrefftzFESpace::GetNDof)
      .def ("SetWavespeed", &TrefftzFESpace::SetWavespeed);
}
#endif // NGS_PYTHON
