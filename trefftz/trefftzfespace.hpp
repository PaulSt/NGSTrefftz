#ifndef FILE_TREFFTZFESPACE_HPP
#define FILE_TREFFTZFESPACE_HPP
#include "trefftzwavefe.hpp"

namespace ngcomp
{

  struct GenericSqrt
  {
    template <typename T> T operator() (T x) const { return sqrt (x); }
    static string Name () { return "sqrt"; }
    void DoArchive (Archive &ar) {}
  };

  class TrefftzFESpace : public FESpace
  {
    int D;
    int fullD;
    int order;
    size_t ndof;
    int nel;
    int nvert;
    int local_ndof;
    float c = 1;
    int useshift = 1;
    int usescale = 1;
    int useqt = 0;
    int heat = 0;
    int heattest = 0;
    int basistype;
    shared_ptr<CoefficientFunction> wavespeedcf = nullptr;
    Matrix<shared_ptr<CoefficientFunction>> GGder;
    Matrix<shared_ptr<CoefficientFunction>> BBder;

  public:
    TrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags);

    void SetWavespeed (shared_ptr<CoefficientFunction> awavespeedcf,
                       shared_ptr<CoefficientFunction> aBBcf = nullptr)
    {
      wavespeedcf = awavespeedcf;
      if (aBBcf || useqt)
        {
          // wavespeedcf = UnaryOpCF(aBBcf/awavespeedcf,GenericSqrt());
          cout << "started auto diff.... ";
          // shared_ptr<CoefficientFunction> GGcf = awavespeedcf;
          // shared_ptr<CoefficientFunction> GGcfx = awavespeedcf;
          shared_ptr<CoefficientFunction> GGcf
              = make_shared<ConstantCoefficientFunction> (1)
                / (awavespeedcf * awavespeedcf);
          shared_ptr<CoefficientFunction> GGcfx
              = make_shared<ConstantCoefficientFunction> (1)
                / (awavespeedcf * awavespeedcf);

          static Timer timerder ("QTrefftzDerivatives");
          static Timer timereval ("QTrefftzDerEval");
          timerder.Start ();
          GGder.SetSize (this->order - 1, (this->order - 2) * (D == 2) + 1);
          for (int ny = 0; ny <= (this->order - 2) * (D == 2); ny++)
            {
              for (int nx = 0; nx <= this->order - 2; nx++)
                {
                  GGder (nx, ny) = GGcfx;
                  GGcfx = GGcfx->Diff (
                      MakeCoordinateCoefficientFunction (0).get (),
                      make_shared<ConstantCoefficientFunction> (1));
                }
              GGcf = GGcf->Diff (MakeCoordinateCoefficientFunction (1).get (),
                                 make_shared<ConstantCoefficientFunction> (1));
              GGcfx = GGcf;
            }
          timerder.Stop ();

          if (!aBBcf)
            {
              aBBcf = make_shared<ConstantCoefficientFunction> (1);
              cout << "SETTING BB TO 1" << endl;
            }
          static Timer timerbb ("QTrefftzBB");
          timerbb.Start ();
          shared_ptr<CoefficientFunction> BBcf = aBBcf;
          shared_ptr<CoefficientFunction> BBcfx = aBBcf;
          BBder.SetSize (this->order, (this->order - 1) * (D == 2) + 1);
          for (int ny = 0; ny <= (this->order - 1) * (D == 2); ny++)
            {
              for (int nx = 0; nx <= this->order - 1; nx++)
                {
                  BBder (nx, ny) = BBcfx;
                  BBcfx = BBcfx->Diff (
                      MakeCoordinateCoefficientFunction (0).get (),
                      make_shared<ConstantCoefficientFunction> (1));
                }
              BBcf = BBcf->Diff (MakeCoordinateCoefficientFunction (1).get (),
                                 make_shared<ConstantCoefficientFunction> (1));
              BBcfx = BBcf;
            }
          timerbb.Stop ();
          cout << "finish" << endl;
        }
    }

    string GetClassName () const override { return "trefftz"; }

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;

    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;

    size_t GetNDof () const override { return ndof; }

    static DocInfo GetDocu ();

  protected:
    template <int D> double Adiam (ElementId ei, double c) const;

    template <int D>
    double Adiam (ElementId ei, shared_ptr<CoefficientFunction> c) const;

    template <int D> Vec<D + 1> ElCenter (ElementId ei) const;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzFESpace (py::module m);
#endif // NGS_PYTHON

#endif
