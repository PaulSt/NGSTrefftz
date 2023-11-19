#include <fem.hpp>
#include <comp.hpp>
#include "scalarmappedfe.hpp"
#include "specialintegrator.hpp"

namespace ngfem
{

  template <int D>
  void GetFacetVertices (ELEMENT_TYPE eltype, int LocalFacetNr,
                         double vertices[][D])
  {
    // int nvert_facet = D==3 ? 4 : D==2 ? 2 : throw Exception("wrong
    // dimension");
    ELEMENT_TYPE etfacet
        = ElementTopology::GetFacetType (eltype, LocalFacetNr);
    int nvert_facet = ElementTopology::GetNVertices (etfacet);
    // int nfaces = ElementTopology::GetNFaces (eltype);
    // int nvert_face = ElementTopology::GetNVertices (faces[LocalFacetNr]);
    auto elvertices = ElementTopology::GetVertices (eltype);
    if (D == 3)
      {
        const FACE *faces = ElementTopology::GetFaces (eltype);
        for (int i = 0; i < nvert_facet; i++)
          {
            for (int j = 0; j < D; j++)
              vertices[i][j] = elvertices[faces[LocalFacetNr][i]][j];
          }
      }
    if (D == 2)
      {
        const EDGE *edges = ElementTopology::GetEdges (eltype);
        for (int i = 0; i < nvert_facet; i++)
          {
            for (int j = 0; j < D; j++)
              vertices[i][j] = elvertices[edges[LocalFacetNr][i]][j];
          }
      }
  }

  template <int D>
  void SpaceTimeDG_FFacetBFI<D>::CalcFacetMatrix (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatMatrix<double> elmat, LocalHeap &lh) const
  {
    static Timer t ("SpaceTimeDG_FFacetBFI ");
    HeapReset hr (lh);
    RegionTracer reg (TaskManager::GetThreadId (), t);

    auto &fel1 = dynamic_cast<const BaseScalarMappedElement &> (volumefel1);
    auto &fel2 = dynamic_cast<const BaseScalarMappedElement &> (volumefel2);
    // int D = volumefel1.Dim ();

    int nd1 = fel1.GetNDof ();
    int nd2 = fel2.GetNDof ();
    if (LocalFacetNr2 == -1)
      nd2 = 0;

    elmat = 0;

    FlatVector<> shape1 (nd1, lh);
    FlatVector<> shape2 (nd2, lh);
    FlatMatrix<> dshape1 (nd1, D, lh);
    FlatMatrix<> dshape2 (nd2, D, lh);
    FlatMatrixFixHeight<2> bmat (nd1 + nd2, lh);
    FlatMatrixFixHeight<2> dbmat (nd1 + nd2, lh);
    Mat<2> dmat;

    ELEMENT_TYPE eltype1 = volumefel1.ElementType ();
    ELEMENT_TYPE eltype2 = volumefel2.ElementType ();
    Facet2ElementTrafo transform1 (eltype1, ElVertices1);
    Facet2ElementTrafo transform2 (eltype2, ElVertices2);
    ELEMENT_TYPE etfacet
        = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    int maxorder = max2 (fel1.Order (), fel2.Order ());
    // IntegrationRule ir (fel1.ElementType (), 2 * maxorder);
    const IntegrationRule &ir_facet
        = SelectIntegrationRule (etfacet, 2 * maxorder);
    ELEMENT_TYPE etffacet = (D == 3 ? ET_SEGM : ET_POINT);
    IntegrationRule ir_ffacet (etffacet, 2 * maxorder);

    const NORMAL *normals1 = ElementTopology::GetNormals (eltype1);
    const NORMAL *normals2 = ElementTopology::GetNormals (eltype2);

    Vec<D> normal_ref1, normal_ref2;
    for (int i = 0; i < D; i++)
      {
        normal_ref1 (i) = normals1[LocalFacetNr1][i];
        normal_ref2 (i) = normals2[LocalFacetNr2][i];
      }

    double vertices[D + 1][D];
    GetFacetVertices<D> (eltype1, LocalFacetNr1, vertices);
    int nvert_facet = ElementTopology::GetNVertices (etfacet);
    double slab_height = numeric_limits<double>::max ();
    double tfix = numeric_limits<double>::max ();
    for (int i = 0; i < nvert_facet; i++)
      {
        FlatVector<> point1 (D, lh);
        FlatVector<> point2 (D, lh);
        eltrans1.CalcPoint (IntegrationPoint (vertices[i], 0), point1);
        eltrans1.CalcPoint (
            IntegrationPoint (vertices[(i + 1) % nvert_facet], 0), point2);
        tfix = min (tfix, point1[D - 1]);
        if (abs (point1[D - 1] - point2[D - 1]) > 0)
          slab_height = min (slab_height, abs (point1[D - 1] - point2[D - 1]));
      }

    for (int l = 0; l < ir_facet.GetNIP () / (maxorder + 1); l++)
      {

        IntegrationPoint ip1 = transform1 (LocalFacetNr1, ir_facet[l]);
        MappedIntegrationPoint<D, D> mip1 (ip1, eltrans1);
        mip1.Point ()[D - 1] = tfix;
        IntegrationPoint ip2 = (LocalFacetNr2 != -1)
                                   ? transform2 (LocalFacetNr2, ir_facet[l])
                                   : ip1;
        MappedIntegrationPoint<D, D> mip2 (ip2, eltrans2);
        mip2.Point ()[D - 1] = tfix;
        // BaseMappedIntegrationPoint *mip;
        // switch (D)
        //{
        // case 2:
        // mip = new MappedIntegrationPoint<2, 2> (ip1, eltrans1);
        // break;
        // case 3:
        // mip = new MappedIntegrationPoint<3, 3> (ip1, eltrans1);
        // break;
        // default:
        // throw Exception ("wrong dimension");
        //}

        Mat<D> inv_jac1 = mip1.GetJacobianInverse ();
        double det1 = mip1.GetJacobiDet ();
        Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;
        double len1 = L2Norm (normal1);
        normal1 /= len1;
        double weight = len1 / slab_height * ir_ffacet[l].Weight ();
        if (abs (normal1[D - 1]) - 1e-8 > 0)
          {
            return;
          }

        fel1.CalcShape (mip1, shape1);
        fel1.CalcDShape (mip1, dshape1);
        double a1 = coef_a->Evaluate (mip1);
        double sig = coef_sig->Evaluate (mip1);
        bmat.Row (0).Range (0, nd1)
            = (LocalFacetNr2 != -1 ? 0.5 : 1.0) * a1 * dshape1 * normal1;
        bmat.Row (1).Range (0, nd1) = shape1;

        if (LocalFacetNr2 != -1)
          {
            Mat<D> inv_jac2 = mip2.GetJacobianInverse ();
            double det2 = mip2.GetJacobiDet ();
            Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;
            double len2 = L2Norm (normal2);
            if (abs (len1 - len2) > 1e-6)
              {
                std::cout << "len :\t" << len1 << "\t=?=\t" << len2
                          << std::endl;
                throw Exception ("DGInnerFacet_LaplaceIntegrator: len1!=len2");
              }
            normal2 /= len2;
            fel2.CalcShape (mip2, shape2);
            fel2.CalcDShape (mip2, dshape2);
            double a2 = coef_a->Evaluate (mip2);
            bmat.Row (0).Range (nd1, nd1 + nd2)
                = -0.5 * a2 * dshape2 * normal2;
            bmat.Row (1).Range (nd1, nd1 + nd2) = -shape2;
          }

        dmat (0, 0) = 0;
        dmat (1, 0) = -1;
        dmat (0, 1) = -1;
        dmat (1, 1) = sig;
        dmat *= weight;
        dbmat = dmat * bmat;
        elmat += Trans (bmat) * dbmat;
      }
  }

  template <int D>
  void SpaceTimeDG_FFacetBFI<D>::CalcFacetMatrix (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &seltrans, FlatArray<int> &SElVertices,
      FlatMatrix<double> elmat, LocalHeap &lh) const
  {
    this->CalcFacetMatrix (volumefel, LocalFacetNr, eltrans, ElVertices,
                           volumefel, -1, eltrans, ElVertices, elmat, lh);
  }

  template <int D>
  void SpaceTimeDG_FFacetBFI<D>::ApplyFacetMatrix (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatVector<double> elx, FlatVector<double> ely, LocalHeap &lh) const
  {
    HeapReset hr (lh);
    auto &fel1 = dynamic_cast<const BaseScalarMappedElement &> (volumefel1);
    auto &fel2 = dynamic_cast<const BaseScalarMappedElement &> (volumefel2);
    int nd1 = fel1.GetNDof ();
    int nd2 = fel2.GetNDof ();
    FlatMatrix<> elmat (nd1 + nd2, nd1 + nd2, lh);
    this->CalcFacetMatrix (volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                           volumefel2, LocalFacetNr2, eltrans2, ElVertices2,
                           elmat, lh);
    ely = elmat * elx;
  }

  template <int D>
  void SpaceTimeDG_FFacetBFI<D>::ApplyFacetMatrix (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &seltrans, FlatArray<int> &SElVertices,
      FlatVector<double> elx, FlatVector<double> ely, LocalHeap &lh) const
  {
    HeapReset hr (lh);
    auto &fel = dynamic_cast<const BaseScalarMappedElement &> (volumefel);
    int nd = fel.GetNDof ();
    FlatMatrix<> elmat (nd, nd, lh);
    this->CalcFacetMatrix (volumefel, LocalFacetNr, eltrans, ElVertices,
                           seltrans, SElVertices, elmat, lh);
    ely = elmat * elx;
  }

  template class SpaceTimeDG_FFacetBFI<2>;
  template class SpaceTimeDG_FFacetBFI<3>;

}

void ExportSpecialIntegrator (py::module m)
{
  using namespace ngfem;
  using namespace ngcomp;
  // py::class_<SpaceTimeDG_FFacetBFI, shared_ptr<SpaceTimeDG_FFacetBFI>,
  // BilinearFormIntegrator> ( m, "SpaceTimeDG_FFacetBFI") .def
  //(py::init<double, shared_ptr<CoefficientFunction>,
  // shared_ptr<CoefficientFunction>> ());
  m.def (
      "SpaceTimeDG_FFacetBFI",
      [] (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> coef_c,
          shared_ptr<CoefficientFunction> coef_sig,
          VorB vb) -> shared_ptr<BilinearFormIntegrator> {
        switch (ma->GetDimension ())
          {
          case 2:
            return make_shared<SpaceTimeDG_FFacetBFI<2>> (coef_c, coef_sig,
                                                          vb);
          case 3:
            return make_shared<SpaceTimeDG_FFacetBFI<3>> (coef_c, coef_sig,
                                                          vb);
          default:
            throw Exception ("wrong dimension");
          }
      },
      py::arg ("mesh"), py::arg ("coef_c"), py::arg ("coef_sig"),
      py::arg ("VOL_or_BND"));

  // py::class_<SpaceTimeDG_FFacetLFI, shared_ptr<SpaceTimeDG_FFacetLFI>,
  // LinearFormIntegrator> ( m, "SpaceTimeDG_FFacetLFI") .def
  //(py::init<shared_ptr<CoefficientFunction>> ());
  // m.def (
  //"SpaceTimeDG_FFacetLFI",
  //[] (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> coef_c,
  // shared_ptr<CoefficientFunction> coef_sig,
  // VorB vb) -> shared_ptr<LinearFormIntegrator> {
  // switch (ma->GetDimension ())
  //{
  // case 2:
  // return make_shared<SpaceTimeDG_FFacetLFI<2>> (ma, coef_c,
  // coef_sig, vb);
  // case 3:
  // return make_shared<SpaceTimeDG_FFacetLFI<3>> (ma, coef_c,
  // coef_sig, vb);
  // default:
  // throw Exception ("wrong dimension");
  //}
  //},
  // py::arg ("mesh"), py::arg ("coef_c"), py::arg ("coef_sig"),
  // py::arg ("VOL_or_BND"));
}
