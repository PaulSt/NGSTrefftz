#include <comp.hpp>
#include "scalarmappedfe.hpp"
#include "specialintegrator.hpp"

namespace ngfem
{

  template <int D>
  bool PointContainedInElement (shared_ptr<ngcomp::MeshAccess> ma,
                                MappedIntegrationPoint<D, D> mip,
                                IntegrationPoint &ip, int elnr)
  {
    netgen::Ngx_Mesh ngmesh = ma->GetNetgenMeshX ();
    // double lamis[3];
    switch (D)
      {
      case 2:
        {
          netgen::Point3d p (mip.GetPoint ()[0], mip.GetPoint ()[1], 0);
          return ngmesh.GetMesh ()->PointContainedIn2DElement (
              p, &ip (0), elnr + 1, false);
        }
      case 3:
        {
          netgen::Point3d p (mip.GetPoint ()[0], mip.GetPoint ()[1],
                             mip.GetPoint ()[2]);
          return ngmesh.GetMesh ()->PointContainedIn3DElement (p, &ip (0),
                                                               elnr + 1);
        }
      default:
        throw Exception ("wrong dimension");
      }
  }

  template <int D>
  int GetElnrNeighbourByPoint (shared_ptr<ngcomp::MeshAccess> ma,
                               MappedIntegrationPoint<D, D> mip,
                               ElementId elnr, IntegrationPoint &neighbourip,
                               int &localfnum)
  {
    auto fnums = ma->GetElFacets (elnr);
    localfnum = 0;
    Array<int> elnums;
    for (int i : fnums)
      {
        ma->GetFacetElements (i, elnums);
        for (size_t elnr2 : elnums)
          {
            if (elnr2 != elnr.Nr ()
                && PointContainedInElement (ma, mip, neighbourip, elnr2))
              return elnr2;
          }
        localfnum++;
      }
    return -1;
  }

  template <int D>
  int GetLocFacetByPoint (shared_ptr<ngcomp::MeshAccess> ma,
                          MappedIntegrationPoint<D, D> mip, int elnr)
  {
    IntegrationPoint ip;
    auto fnums = ma->GetElFacets (elnr);
    for (int i : fnums)
      {
        Array<int> elnums;
        ma->GetFacetElements (i, elnums);
        for (int elnr2 : elnums)
          {
            if (elnr2 != elnr)
              {
                if (PointContainedInElement (ma, mip, ip, elnr2))
                  return elnr2;
              }
          }
      }
  }

  template <int D>
  Vec<D> GetNormal (ELEMENT_TYPE eltyp, int LocalFacetNr,
                    MappedIntegrationPoint<D, D> mip)
  {
    const NORMAL *normals = ElementTopology::GetNormals (eltyp);
    Vec<D> normal_ref;
    for (int i = 0; i < D; i++)
      normal_ref (i) = normals[LocalFacetNr][i];

    Mat<D> inv_jac = mip.GetJacobianInverse ();
    double det = mip.GetJacobiDet ();
    Vec<D> normal = det * Trans (inv_jac) * normal_ref;

    return normal;
  }

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

  void GetSlabInfo (ELEMENT_TYPE eltype, const ElementTransformation &eltrans1,
                    int LocalFacetNr, double &slabheight, double &slabwidth,
                    double &tbot)
  {
    int D = eltrans1.SpaceDim ();
    ELEMENT_TYPE etfacet
        = ElementTopology::GetFacetType (eltype, LocalFacetNr);
    int nvert_facet = ElementTopology::GetNVertices (etfacet);
    // int nfaces = ElementTopology::GetNFaces (eltype);
    // int nvert_face = ElementTopology::GetNVertices (faces[LocalFacetNr]);
    auto elvertices = ElementTopology::GetVertices (eltype);
    const FACE *faces = ElementTopology::GetFaces (eltype);
    const EDGE *edges = ElementTopology::GetEdges (eltype);

    slabheight = numeric_limits<double>::max ();
    slabwidth = D > 2 ? numeric_limits<double>::max () : 1.0;
    tbot = numeric_limits<double>::max ();
    for (int i = 0; i < nvert_facet; i++)
      {
        Vector<> point1 (D);
        if (D == 3)
          {
            IntegrationPoint ip1 (elvertices[faces[LocalFacetNr][i]], 0);
            eltrans1.CalcPoint (ip1, point1);
          }
        else if (D == 2)
          {
            IntegrationPoint ip1 (elvertices[edges[LocalFacetNr][i]], 0);
            eltrans1.CalcPoint (ip1, point1);
          }
        tbot = min (tbot, point1[D - 1]);
        for (int j = 0; j < nvert_facet; j++)
          {
            Vector<> point2 (D);
            if (D == 3)
              {
                IntegrationPoint ip2 (elvertices[faces[LocalFacetNr][j]], 0);
                eltrans1.CalcPoint (ip2, point2);
              }
            else if (D == 2)
              {
                IntegrationPoint ip2 (elvertices[edges[LocalFacetNr][j]], 0);
                eltrans1.CalcPoint (ip2, point2);
              }
            if (abs (point1[D - 1] - point2[D - 1]) > 0)
              slabheight
                  = min (slabheight, abs (point1[D - 1] - point2[D - 1]));
            else if (L2Norm (point1 - point2) > 0)
              slabwidth = min (slabwidth, L2Norm (point1 - point2));
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
        eltrans1.CalcPoint (IntegrationPoint (vertices[i], 0), point1);
        tfix = min (tfix, point1[D - 1]);
        for (int j = 0; j < nvert_facet; j++)
          {
            FlatVector<> point2 (D, lh);
            eltrans1.CalcPoint (IntegrationPoint (vertices[j], 0), point2);
            if (abs (point1[D - 1] - point2[D - 1]) > 0)
              slab_height
                  = min (slab_height, abs (point1[D - 1] - point2[D - 1]));
          }
      }

    for (size_t l = 0; l < ir_facet.GetNIP () / (maxorder + 1); l++)
      {

        IntegrationPoint ip1 = transform1 (LocalFacetNr1, ir_facet[l]);
        MappedIntegrationPoint<D, D> mip1 (ip1, eltrans1);
        mip1.Point ()[D - 1] = tfix;
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
        if (l == 0)
          {
            // cout << "normal1 = " << normal1 << endl;
            // cout << "len1 = " << len1 << endl;
            // cout << "weight = " << weight << endl;
            // cout << "top_or_bottom = " << top_or_bottom << endl;
            // cout << "slab_height = " << slab_height << endl;
            // cout << "tfix = " << tfix << endl;
            // cout << endl;
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
            IntegrationPoint ip2 = transform2 (LocalFacetNr2, ir_facet[l]);
            MappedIntegrationPoint<D, D> mip2 (ip2, eltrans2);
            mip2.Point ()[D - 1] = tfix;
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
      const ElementTransformation &, FlatArray<int> &,
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

  template <int D>
  void SpaceTimeDG_FFacetLFI<D>::CalcFacetVector (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatVector<double> elvec, LocalHeap &lh) const
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

    elvec = 0;

    FlatVector<> shape1 (nd1, lh);
    FlatVector<> shape2 (nd2, lh);
    FlatMatrix<> dshape1 (nd1, D, lh);
    FlatMatrix<> dshape2 (nd2, D, lh);
    FlatMatrixFixHeight<2> bmat (nd1 + nd2, lh);
    Vec<2> bvec = 0;
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
        eltrans1.CalcPoint (IntegrationPoint (vertices[i], 0), point1);
        tfix = min (tfix, point1[D - 1]);
        for (int j = 0; j < nvert_facet; j++)
          {
            FlatVector<> point2 (D, lh);
            eltrans1.CalcPoint (IntegrationPoint (vertices[j], 0), point2);
            if (abs (point1[D - 1] - point2[D - 1]) > 0)
              slab_height
                  = min (slab_height, abs (point1[D - 1] - point2[D - 1]));
          }
      }

    for (size_t l = 0; l < ir_facet.GetNIP () / (maxorder + 1); l++)
      {

        IntegrationPoint ip1 = transform1 (LocalFacetNr1, ir_facet[l]);
        MappedIntegrationPoint<D, D> mip1 (ip1, eltrans1);
        mip1.Point ()[D - 1] = tfix;

        Mat<D> inv_jac1 = mip1.GetJacobianInverse ();
        double det1 = mip1.GetJacobiDet ();
        Vec<D> normal1, normal2;
        normal1 = det1 * Trans (inv_jac1) * normal_ref1;
        double len1 = L2Norm (normal1);
        normal1 /= len1;
        double weight = len1 / slab_height * ir_ffacet[l].Weight ();
        if (abs (normal1[D - 1]) - 1e-8 > 0)
          {
            return;
          }
        if (l == 0)
          {
            // cout << "normal1 = " << normal1 << endl;
            // cout << "len1 = " << len1 << endl;
            // cout << "weight = " << weight << endl;
            // cout << "top_or_bottom = " << top_or_bottom << endl;
            // cout << "slab_height = " << slab_height << endl;
            // cout << "tfix = " << tfix << endl;
            // cout << endl;
          }

        fel1.CalcShape (mip1, shape1);
        fel1.CalcDShape (mip1, dshape1);
        double a1, a2;
        a1 = coef_a->Evaluate (mip1);
        double sig = coef_sig->Evaluate (mip1);

        bmat.Row (0).Range (0, nd1)
            = (LocalFacetNr2 != -1 ? 0.5 : 1.0) * a1 * dshape1 * normal1;
        bmat.Row (1).Range (0, nd1) = shape1;

        if (LocalFacetNr2 != -1)
          {
            IntegrationPoint ip2 = transform2 (LocalFacetNr2, ir_facet[l]);
            MappedIntegrationPoint<D, D> mip2 (ip2, eltrans2);
            mip2.Point ()[D - 1] = tfix;
            Mat<D> inv_jac2 = mip2.GetJacobianInverse ();
            double det2 = mip2.GetJacobiDet ();
            normal2 = det2 * Trans (inv_jac2) * normal_ref2;
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
            a2 = coef_a->Evaluate (mip2);
            bmat.Row (0).Range (nd1, nd1 + nd2)
                = -0.5 * a2 * dshape2 * normal2;
            bmat.Row (1).Range (nd1, nd1 + nd2) = -shape2;
          }

        ////// bvec
        if (LocalFacetNr2 != -1)
          {
            auto gfuh2 = dynamic_pointer_cast<ngcomp::GridFunction> (gfuh);
            shared_ptr<ngcomp::MeshAccess> uhma;
            if (gfuh2)
              {
                uhma = dynamic_pointer_cast<ngcomp::GridFunction> (gfuh)
                           ->GetMeshAccess ();
                // if (!eltrans1.BelongsToMesh (uhma.get ()))
                //{
                IntegrationPoint rip1;
                auto uh_el1 = uhma->FindElementOfPoint (
                    mip1.GetPoint (), rip1,
                    true); // buildtree not yet threadsafe (maybe now ?)
                const ElementTransformation &uh_trafo1
                    = uhma->GetTrafo (uh_el1, lh);
                // BaseMappedIntegrationPoint &uh_mip1 = uh_trafo1 (rip1, lh);
                MappedIntegrationPoint<D, D> uh_mip1 (rip1, uh_trafo1);
                //// ASSUME CONST OUTWARD NORMAL, OR TODO: GET LOCALFACET
                /// NUMBER FOR
                /// LocalFacetNr2 == -1

                IntegrationPoint rip2;
                int uh_LocalFacetNr1, uh_LocalFacetNr2;
                int uh_elnr2 = GetElnrNeighbourByPoint (
                    uhma, mip1, uh_el1, rip2, uh_LocalFacetNr1);
                ElementId uh_el2 (VOL, uh_elnr2);
                IntegrationPoint dip;
                GetElnrNeighbourByPoint (uhma, mip1, uh_el2, dip,
                                         uh_LocalFacetNr2);
                const ElementTransformation &uh_trafo2
                    = uhma->GetTrafo (uh_el2, lh);
                MappedIntegrationPoint<D, D> uh_mip2 (rip2, uh_trafo2);

                Vec<D> uh_normal1 = GetNormal (uhma->GetElType (uh_el1),
                                               uh_LocalFacetNr1, uh_mip1);
                uh_normal1 *= 1 / L2Norm (uh_normal1);

                Vec<D> uh_normal2 = GetNormal (uhma->GetElType (uh_el2),
                                               uh_LocalFacetNr2, uh_mip2);
                uh_normal2 *= 1 / L2Norm (uh_normal2);

                // cout << "uh_mip1: " << uh_mip1.Point () << endl
                //<< "uh_mip2: " << uh_mip2.Point () << endl;
                // cout << "normal1: " << normal1 << endl
                //<< "uh_normal1: " << uh_normal1 << endl
                //<< "normal2: " << normal2 << endl
                //<< "uh_normal2: " << uh_normal2 << endl;

                double uh1 = gfuh->Evaluate (uh_mip1);
                // uh1 = gfuh->Evaluate (uh_trafo1(rip1, lh));
                Vec<D - 1> duh1;
                gfduh->Evaluate (uh_mip1, duh1);

                double uh2 = gfuh->Evaluate (uh_mip2);
                Vec<D - 1> duh2;
                gfduh->Evaluate (uh_mip2, duh2);

                // double test = 0;
                // cout << "uh_normal1" << uh_normal1 << endl;
                // for (int i = 0; i < D - 1; i++)
                // test += 0.5 * (a1*duh1[i] + a2*duh2[i]) * uh_normal1[i] *
                // (i==1);

                // cout << "duh1: " << duh1 << " duh2: " << duh2 << endl;
                // cout << "avg duh n " << test << endl;

                for (int i = 0; i < D - 1; i++)
                  bvec (0)
                      += 0.5 * (a1 * duh1[i] + a2 * duh2[i]) * uh_normal1[i];
                bvec (1) += uh1 - uh2;
                // cout << "gfuh is a gridfunction" << endl;
              }
            else
              {
                IntegrationPoint ip2 = transform2 (LocalFacetNr2, ir_facet[l]);
                MappedIntegrationPoint<D, D> mip2 (ip2, eltrans2);
                double uh1 = gfuh->Evaluate (mip1);
                // uh1 = gfuh->Evaluate (uh_trafo1(rip1, lh));
                Vec<D - 1> duh1;
                gfduh->Evaluate (mip1, duh1);

                double uh2 = gfuh->Evaluate (mip2);
                Vec<D - 1> duh2;
                gfduh->Evaluate (mip2, duh2);

                // double test = 0;
                // cout << "uh_normal1" << uh_normal1 << endl;
                // for (int i = 0; i < D - 1; i++)
                // test += 0.5 * (a1*duh1[i] + a2*duh2[i]) * uh_normal1[i] *
                // (i==1);

                // cout << "duh1: " << duh1 << " duh2: " << duh2 << endl;
                // cout << "avg duh n " << test << endl;

                for (int i = 0; i < D - 1; i++)
                  bvec (0) += 0.5 * (a1 * duh1[i] + a2 * duh2[i]) * normal1[i];
                bvec (1) += uh1 - uh2;
              }
          }
        else
          {
            Vec<D> uh_normal1 = normal1;

            double uh1 = gfuh->Evaluate (mip1);
            Vec<D - 1> duh1;
            gfduh->Evaluate (mip1, duh1);

            for (int i = 0; i < D - 1; i++)
              bvec (0) += (LocalFacetNr2 != -1 ? 0.5 : 1.0) * a1 * duh1[i]
                          * uh_normal1[i];
            bvec (1) += uh1;
          }
        //}

        // cout << "uh1: " << uh << endl;
        //  cout << "duh1: " << duh << endl;

        // cout << "bvec: " << bvec << endl;

        dmat (0, 0) = 0;
        dmat (1, 0) = -1;
        dmat (0, 1) = -1;
        dmat (1, 1) = sig;
        dmat *= weight;
        Vec<2> dbvec = dmat * bvec;
        elvec += Trans (bmat) * dbvec;
      }
    // cout << elvec << endl;
  }

  template <int D>
  void SpaceTimeDG_FFacetLFI<D>::CalcFacetVector (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &, FlatVector<double> elvec,
      LocalHeap &lh) const
  {
    this->CalcFacetVector (volumefel, LocalFacetNr, eltrans, ElVertices,
                           volumefel, -1, eltrans, ElVertices, elvec, lh);
  }

  template class SpaceTimeDG_FFacetLFI<2>;
  template class SpaceTimeDG_FFacetLFI<3>;

  //////////////////////////////////////////////////////////////////////////////
  /// SymbolicFFacet
  //////////////////////////////////////////////////////////////////////////////

  SymbolicFFacetBilinearFormIntegrator ::SymbolicFFacetBilinearFormIntegrator (
      shared_ptr<CoefficientFunction> acf, VorB avb, bool eb)
      : cf (acf), vb (avb), element_boundary (eb)
  {
    simd_evaluate = false;

    if (cf->Dimension () != 1)
      throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append (0);
    test_cum.Append (0);
    cf->TraverseTree ([&] (CoefficientFunction &nodecf) {
      auto proxy = dynamic_cast<ProxyFunction *> (&nodecf);
      if (proxy)
        {
          if (proxy->IsTestFunction ())
            {
              if (!test_proxies.Contains (proxy))
                {
                  test_proxies.Append (proxy);
                  test_cum.Append (test_cum.Last () + proxy->Dimension ());
                }
            }
          else
            {
              if (!trial_proxies.Contains (proxy))
                {
                  trial_proxies.Append (proxy);
                  trial_cum.Append (trial_cum.Last () + proxy->Dimension ());
                }
            }
        }
      else if (nodecf.StoreUserData () && !gridfunction_cfs.Contains (&nodecf))
        gridfunction_cfs.Append (&nodecf);
    });

    neighbor_testfunction = false;
    for (auto proxy : test_proxies)
      if (proxy->IsOther ())
        neighbor_testfunction = true;

    cache_cfs = FindCacheCF (*cf);

    cout << IM (6) << "num test_proxies " << test_proxies.Size () << endl;
    cout << IM (6) << "num trial_proxies " << trial_proxies.Size () << endl;
    cout << IM (6) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM (6) << "cumulated trial_proxy dims " << trial_cum << endl;
  }

  void SymbolicFFacetBilinearFormIntegrator::CalcFacetMatrix (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatMatrix<double> elmat, LocalHeap &lh) const
  {
    T_CalcFacetMatrix (volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                       volumefel2, LocalFacetNr2, eltrans2, ElVertices2, elmat,
                       lh);
  }

  void SymbolicFFacetBilinearFormIntegrator::CalcFacetMatrix (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &seltrans, FlatArray<int> &SElVertices,
      FlatMatrix<double> elmat, LocalHeap &lh) const
  {
    T_CalcFacetMatrix (volumefel, LocalFacetNr, eltrans, ElVertices, seltrans,
                       SElVertices, elmat, lh);
  }

  void SymbolicFFacetBilinearFormIntegrator::CalcFacetMatrix (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatMatrix<Complex> elmat, LocalHeap &lh) const
  {
    if (volumefel1.ComplexShapes () && volumefel2.ComplexShapes ())
      T_CalcFacetMatrix<Complex, Complex> (
          volumefel1, LocalFacetNr1, eltrans1, ElVertices1, volumefel2,
          LocalFacetNr2, eltrans2, ElVertices2, elmat, lh);
    else
      T_CalcFacetMatrix<Complex, double> (
          volumefel1, LocalFacetNr1, eltrans1, ElVertices1, volumefel2,
          LocalFacetNr2, eltrans2, ElVertices2, elmat, lh);
  }

  void SymbolicFFacetBilinearFormIntegrator::CalcFacetMatrix (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &seltrans, FlatArray<int> &SElVertices,
      FlatMatrix<Complex> elmat, LocalHeap &lh) const
  {
    if (volumefel.ComplexShapes ())
      T_CalcFacetMatrix<Complex, Complex> (volumefel, LocalFacetNr, eltrans,
                                           ElVertices, seltrans, SElVertices,
                                           elmat, lh);
    else
      T_CalcFacetMatrix<Complex, double> (volumefel, LocalFacetNr, eltrans,
                                          ElVertices, seltrans, SElVertices,
                                          elmat, lh);
  }

  template <typename TSCAL, typename SCAL_SHAPES>
  void SymbolicFFacetBilinearFormIntegrator ::T_CalcFacetMatrix (
      const FiniteElement &fel1, int LocalFacetNr1,
      const ElementTransformation &trafo1, FlatArray<int> &ElVertices1,
      const FiniteElement &fel2, int LocalFacetNr2,
      const ElementTransformation &trafo2, FlatArray<int> &ElVertices2,
      FlatMatrix<TSCAL> elmat, LocalHeap &lh) const
  {
    elmat = 0.0;

    if (LocalFacetNr2 == -1)
      throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    bool is_mixedfe1 = typeid (fel1) == typeid (const MixedFiniteElement &);
    const MixedFiniteElement *mixedfe1
        = static_cast<const MixedFiniteElement *> (&fel1);
    const FiniteElement &fel1_trial
        = is_mixedfe1 ? mixedfe1->FETrial () : fel1;
    const FiniteElement &fel1_test = is_mixedfe1 ? mixedfe1->FETest () : fel1;
    bool is_mixedfe2 = typeid (fel2) == typeid (const MixedFiniteElement &);
    const MixedFiniteElement *mixedfe2
        = static_cast<const MixedFiniteElement *> (&fel2);
    const FiniteElement &fel2_trial
        = is_mixedfe2 ? mixedfe2->FETrial () : fel2;
    const FiniteElement &fel2_test = is_mixedfe2 ? mixedfe2->FETest () : fel2;

    if (is_mixedfe1 != is_mixedfe2)
      throw Exception ("both sides should have the same type of mixity");

    int maxorder = max2 (max2 (fel1_trial.Order (), fel1_test.Order ()),
                         max2 (fel2_trial.Order (), fel2_test.Order ()));

    auto eltype1 = trafo1.GetElementType ();
    auto eltype2 = trafo2.GetElementType ();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    const IntegrationRule &ir_facet
        = GetIntegrationRule (etfacet, 2 * maxorder + bonus_intorder);
    Facet2ElementTrafo transform1 (eltype1, ElVertices1);
    IntegrationRule &ir_facet_vol1 = transform1 (LocalFacetNr1, ir_facet, lh);

    Facet2ElementTrafo transform2 (eltype2, ElVertices2);
    IntegrationRule &ir_facet_vol2 = transform2 (LocalFacetNr2, ir_facet, lh);

    //// slab info
    double slabheight, slabwidth, tfix;
    GetSlabInfo (eltype1, trafo1, LocalFacetNr1, slabheight, slabwidth, tfix);
    if (slabheight == numeric_limits<double>::max ())
      return;
    int D = trafo1.SpaceDim ();
    ////modify intrule
    ELEMENT_TYPE etffacet = (D == 3 ? ET_SEGM : ET_POINT);
    IntegrationRule ir_ffacet (etffacet, 2 * maxorder);
    size_t ir_size = ir_ffacet.GetNIP ();
    // int ir_size = ir_facet.GetNIP () / (maxorder + 1);
    // cout << "ir_size " << ir_size << " vs " << ir_facet.GetNIP () /
    // (maxorder + 1) << endl;
    ir_facet_vol1.SetSize (ir_size);
    ir_facet_vol2.SetSize (ir_size);
    BaseMappedIntegrationRule &mir1 = trafo1 (ir_facet_vol1, lh);
    BaseMappedIntegrationRule &mir2 = trafo2 (ir_facet_vol2, lh);
    for (size_t i = 0; i < ir_size; i++)
      {
        mir1[i].GetPoint ()[D - 1] = tfix;
        mir2[i].GetPoint ()[D - 1] = tfix;

        double measure1 = slabwidth * ir_ffacet[i].Weight ();
        // mir1[i].GetMeasure () / slabheight // / mir1[i].IP ().Weight ();
        mir1[i].SetMeasure (measure1);
        mir2[i].SetMeasure (measure1);
      }
    // for(int i =0; i < ir_size; i++)
    //{
    // cout << "mir1point " << mir1[i].GetPoint();
    // cout << "mir1[i].GetMeasure() " << mir1[i].GetMeasure() << endl;
    // cout << "mir1[i].GetWeight() " << mir1[i].GetWeight() << endl;
    // cout << "weigt " << ir_ffacet[i].Weight() << endl;
    // cout << "slabheight " << slabheight << " slabwidth " << slabwidth << "
    // facevol " << slabwidth*slabheight << endl; cout << endl << endl;
    //}

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
    mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);

    ProxyUserData ud;
    const_cast<ElementTransformation &> (trafo1).userdata = &ud;

    PrecomputeCacheCF (cache_cfs, mir1, lh);

    for (int k1 : Range (trial_proxies))
      for (int l1 : Range (test_proxies))
        {
          HeapReset hr (lh);

          FlatMatrix<TSCAL> val (mir1.Size (), 1, lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3, TSCAL> proxyvalues (lh, ir_size, proxy2->Dimension (),
                                            proxy1->Dimension ());

          for (size_t k = 0; k < proxy1->Dimension (); k++)
            for (size_t l = 0; l < proxy2->Dimension (); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;

                cf->Evaluate (mir1, val);
                proxyvalues (STAR, l, k) = val.Col (0).Rows (0, ir_size);
              }

          for (size_t i = 0; i < ir_size; i++)
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
            proxyvalues (i, STAR, STAR)
                *= mir1[i].GetMeasure (); // * ir_facet[i].Weight ();

          IntRange trial_range
              = proxy1->IsOther ()
                    ? IntRange (proxy1->Evaluator ()->BlockDim ()
                                    * fel1_trial.GetNDof (),
                                elmat.Width ())
                    : IntRange (0, proxy1->Evaluator ()->BlockDim ()
                                       * fel1_trial.GetNDof ());
          IntRange test_range
              = proxy2->IsOther ()
                    ? IntRange (proxy2->Evaluator ()->BlockDim ()
                                    * fel1_test.GetNDof (),
                                elmat.Height ())
                    : IntRange (0, proxy2->Evaluator ()->BlockDim ()
                                       * fel1_test.GetNDof ());

          auto loc_elmat = elmat.Rows (test_range).Cols (trial_range);
          FlatMatrix<SCAL_SHAPES, ColMajor> bmat1 (proxy1->Dimension (),
                                                   loc_elmat.Width (), lh);
          FlatMatrix<SCAL_SHAPES, ColMajor> bmat2 (proxy2->Dimension (),
                                                   loc_elmat.Height (), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < ir_size; i += BS)
            {
              int rest = min2 (size_t (BS), ir_size - i);
              HeapReset hr (lh);
              FlatMatrix<TSCAL, ColMajor> bdbmat1 (rest * proxy2->Dimension (),
                                                   loc_elmat.Width (), lh);
              FlatMatrix<TSCAL, ColMajor> bbmat2 (rest * proxy2->Dimension (),
                                                  loc_elmat.Height (), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i + j;
                  IntRange r2 = proxy2->Dimension () * IntRange (j, j + 1);
                  if (proxy1->IsOther ())
                    proxy1->Evaluator ()->CalcMatrix (fel2_trial, mir2[ii],
                                                      bmat1, lh);
                  else
                    proxy1->Evaluator ()->CalcMatrix (fel1_trial, mir1[ii],
                                                      bmat1, lh);

                  if (proxy2->IsOther ())
                    proxy2->Evaluator ()->CalcMatrix (fel2_test, mir2[ii],
                                                      bmat2, lh);
                  else
                    proxy2->Evaluator ()->CalcMatrix (fel1_test, mir1[ii],
                                                      bmat2, lh);

                  bdbmat1.Rows (r2) = proxyvalues (ii, STAR, STAR) * bmat1;
                  bbmat2.Rows (r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator ()->UsedDofs (
                  proxy1->IsOther () ? fel2_trial : fel1_trial);
              IntRange r2 = proxy2->Evaluator ()->UsedDofs (
                  proxy2->IsOther () ? fel2_test : fel1_test);
              loc_elmat.Rows (r2).Cols (r1)
                  += Trans (bbmat2.Cols (r2)) * bdbmat1.Cols (r1);
            }
        }
    // cout << "elmat " << elmat << endl;
  }

  template <typename TSCAL, typename SCAL_SHAPES>
  void SymbolicFFacetBilinearFormIntegrator ::T_CalcFacetMatrix (
      const FiniteElement &fel1, int LocalFacetNr1,
      const ElementTransformation &trafo1, FlatArray<int> &ElVertices1,
      const ElementTransformation &strafo, FlatArray<int> &SElVertices1,
      FlatMatrix<TSCAL> elmat, LocalHeap &lh) const
  {
    bool is_mixedfe1 = typeid (fel1) == typeid (const MixedFiniteElement &);
    const MixedFiniteElement *mixedfe1
        = static_cast<const MixedFiniteElement *> (&fel1);
    const FiniteElement &fel1_trial
        = is_mixedfe1 ? mixedfe1->FETrial () : fel1;
    const FiniteElement &fel1_test = is_mixedfe1 ? mixedfe1->FETest () : fel1;

    elmat = 0.0;

    // int maxorder = fel1.Order();
    int maxorder = max2 (fel1_trial.Order (), fel1_test.Order ());

    auto etvol = trafo1.GetElementType ();
    auto etfacet = ElementTopology::GetFacetType (etvol, LocalFacetNr1);

    const IntegrationRule &ir_facet
        = GetIntegrationRule (etfacet, 2 * maxorder + bonus_intorder);
    Facet2ElementTrafo transform1 (etvol, ElVertices1);
    Facet2SurfaceElementTrafo stransform (strafo.GetElementType (),
                                          SElVertices1);

    IntegrationRule &ir_facet_vol1 = transform1 (LocalFacetNr1, ir_facet, lh);
    // IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);  // not yet
    // used ???

    //// slab info
    double slabheight, slabwidth, tfix;
    GetSlabInfo (etvol, trafo1, LocalFacetNr1, slabheight, slabwidth, tfix);
    if (slabheight == numeric_limits<double>::max ())
      return;
    int D = trafo1.SpaceDim ();
    ////modify intrule
    ELEMENT_TYPE etffacet = (D == 3 ? ET_SEGM : ET_POINT);
    IntegrationRule ir_ffacet (etffacet, 2 * maxorder);
    size_t ir_size = ir_ffacet.GetNIP ();
    ir_facet_vol1.SetSize (ir_size);
    BaseMappedIntegrationRule &mir1 = trafo1 (ir_facet_vol1, lh);
    for (size_t i = 0; i < ir_size; i++)
      {
        mir1[i].GetPoint ()[D - 1] = tfix;
        double measure1 = slabwidth * ir_ffacet[i].Weight ();
        mir1[i].SetMeasure (measure1);
      }

    mir1.ComputeNormalsAndMeasure (etvol, LocalFacetNr1);

    ProxyUserData ud;
    const_cast<ElementTransformation &> (trafo1).userdata = &ud;

    PrecomputeCacheCF (cache_cfs, mir1, lh);

    for (int k1 : Range (trial_proxies))
      for (int l1 : Range (test_proxies))
        {
          HeapReset hr (lh);
          FlatMatrix<TSCAL> val (mir1.Size (), 1, lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          if (proxy1->IsOther () || proxy2->IsOther ())
            continue;

          FlatTensor<3, TSCAL> proxyvalues (lh, ir_size, proxy2->Dimension (),
                                            proxy1->Dimension ());

          for (size_t k = 0; k < proxy1->Dimension (); k++)
            for (size_t l = 0; l < proxy2->Dimension (); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;

                cf->Evaluate (mir1, val);
                proxyvalues (STAR, l, k) = val.Col (0).Rows (0, ir_size);
              }

          for (size_t i = 0; i < ir_size; i++)
            proxyvalues (i, STAR, STAR)
                *= mir1[i].GetMeasure (); // * ir_facet[i].Weight ();

          FlatMatrix<SCAL_SHAPES, ColMajor> bmat1 (proxy1->Dimension (),
                                                   elmat.Width (), lh);
          FlatMatrix<SCAL_SHAPES, ColMajor> bmat2 (proxy2->Dimension (),
                                                   elmat.Height (), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < ir_size; i += BS)
            {
              size_t rest = min2 (size_t (BS), ir_size - i);
              HeapReset hr (lh);
              FlatMatrix<TSCAL, ColMajor> bdbmat1 (rest * proxy2->Dimension (),
                                                   elmat.Width (), lh);
              FlatMatrix<TSCAL, ColMajor> bbmat2 (rest * proxy2->Dimension (),
                                                  elmat.Height (), lh);

              for (size_t j = 0; j < rest; j++)
                {
                  size_t ii = i + j;
                  IntRange r2 = proxy2->Dimension () * IntRange (j, j + 1);
                  proxy1->Evaluator ()->CalcMatrix (fel1_trial, mir1[ii],
                                                    bmat1, lh);
                  proxy2->Evaluator ()->CalcMatrix (fel1_test, mir1[ii], bmat2,
                                                    lh);
                  bdbmat1.Rows (r2) = proxyvalues (ii, STAR, STAR) * bmat1;
                  bbmat2.Rows (r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator ()->UsedDofs (fel1_trial);
              IntRange r2 = proxy2->Evaluator ()->UsedDofs (fel1_test);
              elmat.Rows (r2).Cols (r1)
                  += Trans (bbmat2.Cols (r2)) * bdbmat1.Cols (r1) | Lapack;
            }
        }
  }

  SymbolicFFacetLinearFormIntegrator ::SymbolicFFacetLinearFormIntegrator (
      shared_ptr<CoefficientFunction> acf, VorB avb /* , bool eb */)
      : cf (acf), vb (avb) // , element_boundary(eb)
  {
    if (cf->Dimension () != 1)
      throw Exception ("SymblicLFI needs scalar-valued CoefficientFunction");
    test_cum.Append (0);
    cf->TraverseTree ([&] (CoefficientFunction &nodecf) {
      auto proxy = dynamic_cast<ProxyFunction *> (&nodecf);
      if (proxy)
        if (proxy->IsTestFunction ())
          if (!proxies.Contains (proxy))
            {
              proxies.Append (proxy);
              test_cum.Append (test_cum.Last () + proxy->Dimension ());
            }
    });
    cache_cfs = FindCacheCF (*cf);
  }

  template <typename TSCAL>
  void SymbolicFFacetLinearFormIntegrator ::T_CalcFacetVector (
      const FiniteElement &fel1, int LocalFacetNr1,
      const ElementTransformation &trafo1, FlatArray<int> &ElVertices1,
      const FiniteElement &fel2, int LocalFacetNr2,
      const ElementTransformation &trafo2, FlatArray<int> &ElVertices2,
      FlatVector<TSCAL> elvec, LocalHeap &lh) const
  {
    static Timer t ("SymbolicFacetLFI::CalcFacetVector - inner", NoTracing);
    HeapReset hr (lh);

    elvec = 0;

    int maxorder = max2 (fel1.Order (), fel2.Order ());
    auto eltype1 = trafo1.GetElementType ();
    auto eltype2 = trafo2.GetElementType ();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    const IntegrationRule &ir_facet
        = GetIntegrationRule (etfacet, 2 * maxorder + bonus_intorder);
    Facet2ElementTrafo transform1 (eltype1, ElVertices1);
    IntegrationRule &ir_facet_vol1 = transform1 (LocalFacetNr1, ir_facet, lh);
    Facet2ElementTrafo transform2 (eltype2, ElVertices2);
    IntegrationRule &ir_facet_vol2 = transform2 (LocalFacetNr2, ir_facet, lh);
    //// slab info
    double slabheight, slabwidth, tfix;
    GetSlabInfo (eltype1, trafo1, LocalFacetNr1, slabheight, slabwidth, tfix);
    if (slabheight == numeric_limits<double>::max ())
      return;
    int D = trafo1.SpaceDim ();
    ////modify intrule
    ELEMENT_TYPE etffacet = (D == 3 ? ET_SEGM : ET_POINT);
    IntegrationRule ir_ffacet (etffacet, 2 * maxorder);
    size_t ir_size = ir_ffacet.GetNIP ();
    ir_facet_vol1.SetSize (ir_size);
    ir_facet_vol2.SetSize (ir_size);
    BaseMappedIntegrationRule &mir1 = trafo1 (ir_facet_vol1, lh);
    BaseMappedIntegrationRule &mir2 = trafo2 (ir_facet_vol2, lh);
    for (size_t i = 0; i < ir_size; i++)
      {
        mir1[i].GetPoint ()[D - 1] = tfix;
        mir2[i].GetPoint ()[D - 1] = tfix;
        double measure1 = slabwidth * ir_ffacet[i].Weight ();
        mir1[i].SetMeasure (measure1);
        mir2[i].SetMeasure (measure1);
      }

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);

    ProxyUserData ud;
    const_cast<ElementTransformation &> (trafo1).userdata = &ud;

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
    mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);

    PrecomputeCacheCF (cache_cfs, mir1, lh);
    FlatMatrix<TSCAL> val (ir_facet.Size (), 1, lh);
    for (auto proxy : proxies)
      {
        HeapReset hr (lh);
        FlatMatrix<TSCAL> proxyvalues (ir_size, proxy->Dimension (), lh);

        for (size_t k = 0; k < proxy->Dimension (); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf->Evaluate (mir1, val); // needed for grad(u), mesh_size, but:
                                      // index = material index
            proxyvalues.Col (k) = val.Col (0).Rows (0, ir_size);
          }

        for (size_t i = 0; i < ir_size; i++)
          proxyvalues.Row (i)
              *= mir1[i].GetMeasure (); // * ir_facet[i].Weight(); // use also
                                        // smir here ?

        IntRange range = proxy->IsOther ()
                             ? IntRange (proxy->Evaluator ()->BlockDim ()
                                             * fel1.GetNDof (),
                                         elvec.Size ())
                             : IntRange (0, proxy->Evaluator ()->BlockDim ()
                                                * fel1.GetNDof ());

        FlatVector<TSCAL> localvec (range.Size (), lh);
        localvec = 0.0;

        if (proxy->IsOther ())
          proxy->Evaluator ()->ApplyTrans (fel2, mir2, proxyvalues, localvec,
                                           lh);
        else
          proxy->Evaluator ()->ApplyTrans (fel1, mir1, proxyvalues, localvec,
                                           lh);
        elvec.Range (range) += localvec;
      }
  }

  template <typename TSCAL>
  void SymbolicFFacetLinearFormIntegrator ::T_CalcFacetVector (
      const FiniteElement &fel1, int LocalFacetNr,
      const ElementTransformation &trafo1, FlatArray<int> &ElVertices,
      const ElementTransformation &strafo, FlatVector<TSCAL> elvec,
      LocalHeap &lh) const
  {
    static Timer t ("SymbolicFacetLFI::CalcFacetVector - boundary", NoTracing);
    HeapReset hr (lh);

    elvec = 0;

    FlatVector<TSCAL> elvec1 (elvec.Size (), lh);

    int maxorder = fel1.Order ();

    auto eltype1 = trafo1.GetElementType ();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    const IntegrationRule &ir_facet
        = GetIntegrationRule (etfacet, 2 * maxorder + bonus_intorder);
    Facet2ElementTrafo transform1 (eltype1, ElVertices);
    IntegrationRule &ir_facet_vol1 = transform1 (LocalFacetNr, ir_facet, lh);
    auto &smir = strafo (ir_facet, lh);
    //// slab info
    double slabheight, slabwidth, tfix;
    GetSlabInfo (eltype1, trafo1, LocalFacetNr, slabheight, slabwidth, tfix);
    if (slabheight == numeric_limits<double>::max ())
      return;
    int D = trafo1.SpaceDim ();
    ////modify intrule
    ELEMENT_TYPE etffacet = (D == 3 ? ET_SEGM : ET_POINT);
    IntegrationRule ir_ffacet (etffacet, 2 * maxorder);
    size_t ir_size = ir_ffacet.GetNIP ();
    ir_facet_vol1.SetSize (ir_size);
    BaseMappedIntegrationRule &mir1 = trafo1 (ir_facet_vol1, lh);
    for (size_t i = 0; i < ir_size; i++)
      {
        mir1[i].GetPoint ()[D - 1] = tfix;
        double measure1 = slabwidth * ir_ffacet[i].Weight ();
        mir1[i].SetMeasure (measure1);
      }

    mir1.SetOtherMIR (&smir);
    smir.SetOtherMIR (&mir1);

    // evaluate proxy-values
    ProxyUserData ud;
    const_cast<ElementTransformation &> (trafo1).userdata = &ud;
    const_cast<ElementTransformation &> (strafo).userdata = &ud;

    RegionTimer reg (t);

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr);

    PrecomputeCacheCF (cache_cfs, mir1, lh);
    FlatMatrix<TSCAL> val (ir_size, 1, lh);
    for (auto proxy : proxies)
      {
        HeapReset hr (lh);
        FlatMatrix<TSCAL> proxyvalues (ir_size, proxy->Dimension (), lh);

        for (size_t k = 0; k < proxy->Dimension (); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf->Evaluate (mir1, val); // needed for grad(u), mesh_size, but:
                                      // index = material index
            // cf -> Evaluate (smir, val);
            proxyvalues.Col (k) = val.Col (0).Rows (0, ir_size);
          }

        for (size_t i = 0; i < ir_size; i++)
          proxyvalues.Row (i)
              *= mir1[i].GetMeasure (); // * ir_facet[i].Weight(); // use also
                                        // smir here ?

        elvec1 = 0.0;
        proxy->Evaluator ()->ApplyTrans (fel1, mir1, proxyvalues, elvec1, lh);
        elvec += elvec1;
      }
  }

  void SymbolicFFacetLinearFormIntegrator::CalcFacetVector (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatVector<double> elvec, LocalHeap &lh) const
  {
    T_CalcFacetVector (volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                       volumefel2, LocalFacetNr2, eltrans2, ElVertices2, elvec,
                       lh);
  }

  void SymbolicFFacetLinearFormIntegrator::CalcFacetVector (
      const FiniteElement &volumefel1, int LocalFacetNr1,
      const ElementTransformation &eltrans1, FlatArray<int> &ElVertices1,
      const FiniteElement &volumefel2, int LocalFacetNr2,
      const ElementTransformation &eltrans2, FlatArray<int> &ElVertices2,
      FlatVector<Complex> elvec, LocalHeap &lh) const
  {
    T_CalcFacetVector (volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                       volumefel2, LocalFacetNr2, eltrans2, ElVertices2, elvec,
                       lh);
  }

  void SymbolicFFacetLinearFormIntegrator::CalcFacetVector (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &seltrans, FlatVector<double> elvec,
      LocalHeap &lh) const
  {
    T_CalcFacetVector (volumefel, LocalFacetNr, eltrans, ElVertices, seltrans,
                       elvec, lh);
  }

  void SymbolicFFacetLinearFormIntegrator::CalcFacetVector (
      const FiniteElement &volumefel, int LocalFacetNr,
      const ElementTransformation &eltrans, FlatArray<int> &ElVertices,
      const ElementTransformation &seltrans, FlatVector<Complex> elvec,
      LocalHeap &lh) const
  {
    T_CalcFacetVector (volumefel, LocalFacetNr, eltrans, ElVertices, seltrans,
                       elvec, lh);
  }

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
  // LinearFormIntegrator> (m, "SpaceTimeDG_FFacetLFI")
  //.def (py::init<shared_ptr<CoefficientFunction>> ());
  m.def (
      "SpaceTimeDG_FFacetLFI",
      [] (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> gfuh,
          shared_ptr<CoefficientFunction> gfduh,
          shared_ptr<CoefficientFunction> coef_c,
          shared_ptr<CoefficientFunction> coef_sig,
          VorB vb) -> shared_ptr<LinearFormIntegrator> {
        switch (ma->GetDimension ())
          {
          case 2:
            return make_shared<SpaceTimeDG_FFacetLFI<2>> (
                ma, gfuh, gfduh, coef_c, coef_sig, vb);
          case 3:
            return make_shared<SpaceTimeDG_FFacetLFI<3>> (
                ma, gfuh, gfduh, coef_c, coef_sig, vb);
          default:
            throw Exception ("wrong dimension");
          }
      },
      py::arg ("mesh"), py::arg ("gfuh"), py::arg ("gfduh"),
      py::arg ("coef_c"), py::arg ("coef_sig"), py::arg ("VOL_or_BND"));

  m.def (
      "FFacetBFI",
      [] (shared_ptr<CoefficientFunction> cf, VorB vb, bool element_boundary,
          bool skeleton, optional<variant<Region, py::list>> definedon,
          IntegrationRule ir, int bonus_intorder,
          shared_ptr<BitArray> definedonelem, bool simd_evaluate,
          VorB element_vb, bool geom_free,
          shared_ptr<GridFunction> deformation) {
        if (definedon.has_value ())
          if (auto defregion = get_if<Region> (&*definedon); defregion)
            vb = VorB (*defregion);

        if (element_boundary)
          element_vb = BND;
        // check for DG terms
        bool has_other = false;

        cf->TraverseTree ([&has_other] (CoefficientFunction &cf) {
          if (dynamic_cast<ProxyFunction *> (&cf))
            if (dynamic_cast<ProxyFunction &> (cf).IsOther ())
              has_other = true;
        });
        // if (has_other && !element_boundary && !skeleton)
        if (has_other && (element_vb != BND) && !skeleton)
          throw Exception ("DG-facet terms need either skeleton=True or "
                           "element_boundary=True");

        shared_ptr<BilinearFormIntegrator> bfi;
        if (!has_other && !skeleton)
          throw Exception ("DG-ffacet terms need skeleton=True");
        // bfi = make_shared<SymbolicBilinearFormIntegrator> (cf, vb,
        // element_vb);
        else
          bfi = make_shared<SymbolicFFacetBilinearFormIntegrator> (
              cf, vb, element_boundary);
        bfi->geom_free = geom_free;
        if (definedon.has_value ())
          {
            if (auto defpylist = get_if<py::list> (&*definedon); defpylist)
              {
                Array<int> defon = makeCArray<int> (*defpylist);
                for (int &d : defon)
                  d--;
                bfi->SetDefinedOn (defon);
              }
            if (auto defregion = get_if<Region> (&*definedon); defregion)
              bfi->SetDefinedOn (defregion->Mask ());
          }
        bfi->SetBonusIntegrationOrder (bonus_intorder);
        if (ir.Size ())
          {
            cout << IM (1)
                 << "WARNING: Setting the integration rule for all element "
                    "types is deprecated, use "
                    "BFI.SetIntegrationRule(ELEMENT_TYPE, IntegrationRule) "
                    "instead!"
                 << endl;
            /*
            dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (bfi)
              ->SetIntegrationRule(ir);
            */
            bfi->SetIntegrationRule (ir);
          }

        bfi->SetSimdEvaluate (simd_evaluate);
        bfi->SetDeformation (deformation);
        if (definedonelem)
          bfi->SetDefinedOnElements (definedonelem);
        return shared_ptr<BilinearFormIntegrator> (bfi);
      },
      py::arg ("form"), py::arg ("VOL_or_BND") = VOL,
      py::arg ("element_boundary") = false, py::arg ("skeleton") = true,
      py::arg ("definedon") = nullptr,
      py::arg ("intrule") = IntegrationRule (), py::arg ("bonus_intorder") = 0,
      py::arg ("definedonelements") = nullptr,
      py::arg ("simd_evaluate") = false, py::arg ("element_vb") = VOL,
      py::arg ("geom_free") = false,
      py::arg ("deformation") = shared_ptr<GridFunction> (),
      docu_string (R"raw_string(
A symbolic bilinear form integrator, operating on facet intersections.

Parameters:

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND

skeleton : bool
  must be True
)raw_string"));

  m.def (
      "FFacetLFI",
      [] (shared_ptr<CoefficientFunction> cf, VorB vb, bool, bool skeleton,
          optional<variant<Region, py::list>> definedon, IntegrationRule ir,
          int bonus_intorder, shared_ptr<BitArray> definedonelem,
          bool simd_evaluate, VorB, shared_ptr<GridFunction> deformation) {
        if (definedon.has_value ())
          if (auto defregion = get_if<Region> (&*definedon); defregion)
            vb = VorB (*defregion);

        shared_ptr<LinearFormIntegrator> lfi;
        if (!skeleton)
          // lfi = make_shared<SymbolicLinearFormIntegrator> (cf, vb,
          // element_vb);
          throw Exception ("DG-ffacet terms need skeleton=True");
        else
          lfi = make_shared<SymbolicFFacetLinearFormIntegrator> (
              cf, vb /* , element_boundary */);

        if (definedon.has_value ())
          {
            if (auto defpylist = get_if<py::list> (&*definedon); defpylist)
              {
                Array<int> defon = makeCArray<int> (*defpylist);
                for (int &d : defon)
                  d--;
                lfi->SetDefinedOn (defon);
              }
            if (auto defregion = get_if<Region> (&*definedon); defregion)
              lfi->SetDefinedOn (defregion->Mask ());
          }

        lfi->SetSimdEvaluate (simd_evaluate);
        lfi->SetDeformation (deformation);

        lfi->SetBonusIntegrationOrder (bonus_intorder);
        if (ir.Size ())
          {
            cout << IM (1)
                 << "WARNING: Setting the integration rule for all element "
                    "types is deprecated, use "
                    "LFI.SetIntegrationRule(ELEMENT_TYPE, IntegrationRule) "
                    "instead!"
                 << endl;
            dynamic_pointer_cast<SymbolicLinearFormIntegrator> (lfi)
                ->SetIntegrationRule (ir);
          }

        if (definedonelem)
          lfi->SetDefinedOnElements (definedonelem);

        return shared_ptr<LinearFormIntegrator> (lfi);
      },
      py::arg ("form"), py::arg ("VOL_or_BND") = VOL,
      py::arg ("element_boundary") = false, py::arg ("skeleton") = true,
      py::arg ("definedon") = nullptr,
      py::arg ("intrule") = IntegrationRule (), py::arg ("bonus_intorder") = 0,
      py::arg ("definedonelements") = nullptr,
      py::arg ("simd_evaluate") = false, py::arg ("element_vb") = VOL,
      py::arg ("deformation") = shared_ptr<GridFunction> (),
      docu_string (R"raw_string(
A symbolic linear form integrator, operating on facet intersections.

Parameters:

VOL_or_BND : ngsolve.comp.VorB
  input VOL, BND

skeleton : bool
  must be True
)raw_string"));
}
