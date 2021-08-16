#ifndef SPECIALCOEFFICIENTFUNCTION_HPP
#define SPECIALCOEFFICIENTFUNCTION_HPP

#include <fem.hpp>
#include <comp.hpp>
#include <multigrid.hpp>
#include <h1lofe.hpp>
#include <regex>

using namespace ngcomp;
namespace ngfem
{

  class ClipCoefficientFunction : public CoefficientFunction
  {
  private:
    shared_ptr<CoefficientFunction> coef;
    double clipvalue;
    int clipdim;

  public:
    ///
    ClipCoefficientFunction (shared_ptr<CoefficientFunction> acoef,
                             int adimension, int aclipdim, double aclipvalue,
                             bool ais_complex = false)
        : CoefficientFunction (adimension, ais_complex), coef (acoef),
          clipdim (aclipdim), clipvalue (aclipvalue)
    {
      ;
    }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint &ip) const;
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule &ir,
                           FlatMatrix<double> values) const;
    virtual void EvaluateStdRule (const BaseMappedIntegrationRule &ir,
                                  FlatMatrix<double> values) const;
  };

  class IntegrationPointFunction : public CoefficientFunction
  {
  public:
    IntegrationPointFunction (shared_ptr<MeshAccess> mesh,
                              IntegrationRule &intrule, Vector<> ipdata)
        : CoefficientFunction (1)
    {
      values.resize (mesh->GetNE ());
      int elnr = 0;
      for (auto &vec : values)
        {
          vec.resize (intrule.GetNIP ());
          for (int i = 0; i < vec.size (); i++)
            {
              // input data from vector with mip values sorted per element
              vec[i] = ipdata[intrule.Size () * elnr + i];
            }
          elnr++;
        }
    }

    IntegrationPointFunction (shared_ptr<MeshAccess> mesh,
                              IntegrationRule &intrule, Matrix<> ipdata)
        : CoefficientFunction (1)
    {
      values.resize (mesh->GetNE ());
      int elnr = 0;
      for (auto &vec : values)
        {
          vec.resize (intrule.GetNIP ());
          for (int i = 0; i < vec.size (); i++)
            {
              // input data from matrix with elnr in rows, mip values in cols
              vec[i] = ipdata (elnr, i);
            }
          elnr++;
        }
    }

    virtual double Evaluate (const BaseMappedIntegrationPoint &ip) const
    {
      int p = ip.GetIPNr ();
      int el = ip.GetTransformation ().GetElementNr ();

      if (p < 0 || p >= values[el].size ())
        {
          cout << "got illegal integration point number " << p << endl;
          return 0;
        }

      return values[el][p];
    }

    void PrintTable ()
    {
      for (int i = 0; i < values.size (); i++)
        {
          for (int j = 0; j < values[i].size (); j++)
            {
              cout << values[i][j] << ", ";
            }
          cout << endl;
        }
      cout << endl;
    }

  private:
    vector<vector<double>> values;
  };

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

  class WeightedRadiusFunction : public CoefficientFunction
  {
  private:
    Vector<> values;

  public:
    WeightedRadiusFunction (shared_ptr<MeshAccess> mesh,
                            shared_ptr<CoefficientFunction> wavespeedcf)
        : CoefficientFunction (1)
    {
      LocalHeap lh (1000 * 1000 * 100);
      values.SetSize (mesh->GetNE ());
      int elnr = 0;
      for (auto &vec : values)
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
              double c1, c2;
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
          elnr++;
        }
    }

    virtual double Evaluate (const BaseMappedIntegrationPoint &ip) const
    {
      int p = ip.GetIPNr ();
      int el = ip.GetTransformation ().GetElementNr ();

      if (el < 0 || el >= values.Size ())
        {
          cout << "got illegal element number " << el << endl;
          return 0;
        }

      return values[el];
    }
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSpecialCoefficientFunction (py::module m);
#endif // NGS_PYTHON

#endif
