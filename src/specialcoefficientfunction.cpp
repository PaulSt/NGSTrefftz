#include "specialcoefficientfunction.hpp"
#include <fem.hpp>
#include <comp.hpp>
#include <h1lofe.hpp>
#include <multigrid.hpp>

namespace ngfem
{
  using namespace ngcomp;

  double ClipCoefficientFunction ::Evaluate (
      const BaseMappedIntegrationPoint &ip) const
  {
    ip.GetPoint ()[clipdim] = clipvalue;
    return coef->Evaluate (ip);
  }

  void
  ClipCoefficientFunction ::Evaluate (const BaseMappedIntegrationPoint &ip,
                                      FlatVector<> result) const
  {
    ip.GetPoint ()[clipdim] = clipvalue;
    coef->Evaluate (ip, result);
  }

  IntegrationPointFunction ::IntegrationPointFunction (
      shared_ptr<MeshAccess> mesh, IntegrationRule &intrule, Vector<> ipdata)
      : CoefficientFunction (1)
  {
    this->ma = mesh;
    this->intrule.SetSize (intrule.Size ());
    for (size_t i = 0; i < intrule.Size (); i++)
      this->intrule[i] = (intrule)[i];
    this->intrule.SetDim (intrule.Dim ());

    values.resize (mesh->GetNE ());
    int elnr = 0;
    for (auto &vec : values)
      {
        vec.resize (intrule.GetNIP ());
        for (size_t i = 0; i < vec.size (); i++)
          {
            // input data from vector with mip values sorted per element
            vec[i] = ipdata[intrule.Size () * elnr + i];
          }
        elnr++;
      }
  }

  IntegrationPointFunction ::IntegrationPointFunction (
      shared_ptr<MeshAccess> mesh, IntegrationRule &intrule, Matrix<> ipdata)
      : CoefficientFunction (1)
  {
    this->ma = mesh;
    this->intrule.SetSize (intrule.Size ());
    for (size_t i = 0; i < intrule.Size (); i++)
      this->intrule[i] = (intrule)[i];
    this->intrule.SetDim (intrule.Dim ());

    values.resize (mesh->GetNE ());
    int elnr = 0;
    for (auto &vec : values)
      {
        vec.resize (intrule.GetNIP ());
        for (size_t i = 0; i < vec.size (); i++)
          {
            // input data from matrix with elnr in rows, mip values in cols
            vec[i] = ipdata (elnr, i);
          }
        elnr++;
      }
  }

  double IntegrationPointFunction ::Evaluate (
      const BaseMappedIntegrationPoint &ip) const
  {
    size_t p = ip.GetIPNr ();
    int el = ip.GetTransformation ().GetElementNr ();

    if (p >= values[el].size ())
      {
        cout << "got illegal integration point number " << p << endl;
        return 0;
      }

    return values[el][p];
  }

  void IntegrationPointFunction ::PrintTable ()
  {
    for (size_t i = 0; i < values.size (); i++)
      {
        for (size_t j = 0; j < values[i].size (); j++)
          {
            cout << values[i][j] << ", ";
          }
        cout << endl;
      }
    cout << endl;
  }

  vector<vector<double>> IntegrationPointFunction ::Export ()
  {
    LocalHeap lh (1000 * 1000 * 100, "export intpointfct");

    vector<vector<double>> pointnvalues;
    pointnvalues.resize (ma->GetNE () * intrule.GetNIP ());
    for (size_t elnr = 0; elnr < ma->GetNE (); elnr++)
      {
        switch (ma->GetDimension ())
          {
          case 2:
            {
              MappedIntegrationRule<2, 2> mir (intrule,
                                               ma->GetTrafo (elnr, lh), lh);
              for (auto mip : mir)
                {
                  auto &vec = pointnvalues[elnr * intrule.GetNIP ()
                                           + mip.GetIPNr ()];
                  vec.resize (3);
                  for (int i = 0; i < 2; i++)
                    vec[i] = mip.Point ()[i];
                  vec[2] = values[elnr][mip.GetIPNr ()];
                }
              break;
            }
          case 3:
            {
              MappedIntegrationRule<3, 3> mir (intrule,
                                               ma->GetTrafo (elnr, lh), lh);
              for (auto mip : mir)
                {
                  auto &vec = pointnvalues[elnr * intrule.GetNIP ()
                                           + mip.GetIPNr ()];
                  vec.resize (4);
                  for (int i = 0; i < 3; i++)
                    vec[i] = mip.Point ()[i];
                  vec[3] = values[elnr][mip.GetIPNr ()];
                }
              break;
            }
          }
      }
    return pointnvalues;
  }

  WeightedRadiusFunction ::WeightedRadiusFunction (
      shared_ptr<MeshAccess> mesh, shared_ptr<CoefficientFunction> wavespeedcf)
      : CoefficientFunction (1)
  {
    LocalHeap lh (1000 * 1000 * 100);
    values.SetSize (mesh->GetNE ());
    for (size_t elnr = 0; elnr < mesh->GetNE (); elnr++)
      {
        double anisotropicdiam = 0.0;
        int D = mesh->GetDimension ();
        Vector<> center (D);
        center = 0;
        Vector<> v1 (D);
        auto vertices_index = mesh->GetElVertices (elnr);
        switch (mesh->GetDimension ())
          {
          case 2:
            for (auto vertex : vertices_index)
              center += mesh->GetPoint<2> (vertex);
            break;
          case 3:
            for (auto vertex : vertices_index)
              center += mesh->GetPoint<3> (vertex);
            break;
          }
        center *= (1.0 / vertices_index.Size ());

        ElementId ei = ElementId (elnr);
        IntegrationRule ir (mesh->GetElType (ei), 0);
        ElementTransformation &trafo = mesh->GetTrafo (ei, lh);

        double maxc = 0;
        for (auto vertex1 : vertices_index)
          {
            double c1 = 0, c2 = 0;
            switch (mesh->GetDimension ())
              {
              case 2:
                {
                  MappedIntegrationPoint<2, 2> mip (ir[0], trafo);
                  v1 = mesh->GetPoint<2> (vertex1);
                  mip.Point () = v1;
                  c1 = wavespeedcf->Evaluate (mip);
                  maxc = max (abs (c1), maxc);
                  mip.Point () = center;
                  c2 = wavespeedcf->Evaluate (mip);
                  maxc = max (abs (c2), maxc);
                  break;
                }
              case 3:
                {
                  MappedIntegrationPoint<3, 3> mip (ir[0], trafo);
                  v1 = mesh->GetPoint<3> (vertex1);
                  mip.Point () = v1;
                  c1 = wavespeedcf->Evaluate (mip);
                  maxc = max (abs (c1), maxc);
                  mip.Point () = center;
                  c2 = wavespeedcf->Evaluate (mip);
                  maxc = max (abs (c2), maxc);
                  break;
                }
              }
            anisotropicdiam = max (
                anisotropicdiam,
                sqrt (L2Norm2 (v1 (0, D - 1) - center (0, D - 1))
                      + pow (c1 * v1 (D - 1) - c2 * center (D - 1), 2)));
          }
        values[elnr] = anisotropicdiam / maxc;
      }
  }

  double WeightedRadiusFunction ::Evaluate (
      const BaseMappedIntegrationPoint &ip) const
  {
    size_t el = ip.GetTransformation ().GetElementNr ();

    if (el >= values.Size ())
      {
        cout << "got illegal element number " << el << endl;
        return 0;
      }

    return values[el];
  }

  // class TrefftzCoefficientFunction : public CoefficientFunction
  //{
  // int basisfunction;
  // TrefftzWaveFE<3> treff = TrefftzWaveFE<3>(4,1,ET_TRIG,0);

  // public:
  // TrefftzCoefficientFunction()
  //: CoefficientFunction(1) { ; }

  // TrefftzCoefficientFunction(int basis)
  //: CoefficientFunction(1) { basisfunction = basis; }

  // virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const
  // override
  //{
  // FlatVector<double> point = mip.GetPoint();

  // int ndof = treff.GetNDof();
  // cout  << "nr: " << basisfunction << " / " << ndof << endl;
  // Vector<double> shape(ndof);
  ////Matrix<double> shape(ndof,2);
  // treff.CalcShape(mip,shape);
  // return shape[basisfunction];
  //}
  //};

  double PrintCF::Evaluate (const BaseMappedIntegrationPoint &mip) const
  {
    if (ofs)
      {
        for (int i = 0; i < mip.GetTransformation ().SpaceDim (); i++)
          {
            if (i > 0)
              *ofs << "\t";
            *ofs << mip.GetPoint ()[i];
          }
        *ofs << "\t" << mip.GetWeight ();
        *ofs << endl;
      }
    return 1.0;
  }

}

#ifdef NGS_PYTHON
void ExportSpecialCoefficientFunction (py::module m)
{
  using namespace ngcomp;
  typedef shared_ptr<CoefficientFunction> spCF;

  m.def (
      "ClipCoefficientFunction",
      [] (spCF cf_x, int aclipdim,
          double aclipvalue) -> shared_ptr<CoefficientFunction> {
        auto pcf = make_shared<ClipCoefficientFunction> (
            cf_x, cf_x->Dimension (), aclipdim, aclipvalue, false);
        pcf->SetDimension (pcf->Dimension ());
        return pcf;
      },
      py::call_guard<py::gil_scoped_release> ());

  py::class_<IntegrationPointFunction, shared_ptr<IntegrationPointFunction>,
             CoefficientFunction> (m, "IntegrationPointFunction")
      .def (py::init ([] (shared_ptr<MeshAccess> mesh,
                          IntegrationRule &intrule, Vector<> data) {
              return new IntegrationPointFunction (mesh, intrule, data);
            }),
            py::arg ("mesh"), py::arg ("intrule"), py::arg ("Vector"))
      .def (py::init ([] (shared_ptr<MeshAccess> mesh,
                          IntegrationRule &intrule, Matrix<> data) {
              return new IntegrationPointFunction (mesh, intrule, data);
            }),
            py::arg ("mesh"), py::arg ("intrule"), py::arg ("Matrix"))
      .def ("PrintTable", &IntegrationPointFunction::PrintTable)
      .def ("Export", &IntegrationPointFunction::Export);

  py::class_<WeightedRadiusFunction, shared_ptr<WeightedRadiusFunction>,
             CoefficientFunction> (m, "WeightedRadiusFunction")
      .def (py::init ([] (shared_ptr<MeshAccess> mesh,
                          shared_ptr<CoefficientFunction> wavespeedcf) {
              return new WeightedRadiusFunction (mesh, wavespeedcf);
            }),
            py::arg ("mesh"), py::arg ("CoefficientFunction"));

  py::class_<AdjacentFaceSizeCF, shared_ptr<AdjacentFaceSizeCF>,
             CoefficientFunction> (m, "AdjacentFaceSizeCF")
      .def (py::init ([] () { return new AdjacentFaceSizeCF (); }));

  // py::class_<TrefftzCoefficientFunction,
  // shared_ptr<TrefftzCoefficientFunction>, CoefficientFunction> (m,
  //"TrefftzCoefficient", "") .def(py::init<>()) .def(py::init<int>())
  //;

  py::class_<PrintCF, shared_ptr<PrintCF>, CoefficientFunction> (
      m, "PrintCF", docu_string (R"raw_string(
CoefficientFunction that writes integration point (in world coords.)
into a file whenever evaluated at one.
)raw_string"))
      .def (py::init ([] (const string &a_filename) -> shared_ptr<PrintCF> {
              return make_shared<PrintCF> (a_filename);
            }),
            py::arg ("filename"), docu_string (R"raw_string(
Constructor of PrintCF (integration point (in world coords.) printing coefficientfunction).
  Argument: filename (string) : name of the file where the values shall be printed
)raw_string"));
}
#endif // NGS_PYTHON
