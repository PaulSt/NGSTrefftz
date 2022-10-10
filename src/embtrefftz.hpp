#ifndef FILE_SVDTREFFTZ_HPP
#define FILE_SVDTREFFTZ_HPP
#include <comp.hpp>
#include <python_comp.hpp>
#include <fem.hpp>
#include <integratorcf.hpp>
#include <variant>
#include <bla.hpp>

namespace ngcomp
{
  template <class SCAL>
  std::tuple<shared_ptr<BaseMatrix>, shared_ptr<BaseVector>>
  EmbTrefftz (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
              shared_ptr<SumOfIntegrals> lf, double eps,
              shared_ptr<FESpace> fes_test, int tndof,
              std::map<std::string, Vector<SCAL>> *stats = nullptr);

  class EmbTrefftzFESpace
      : public L2HighOrderFESpace //, public
                                  //std::enable_shared_from_this<EmbTrefftzFESpace>
  {
    shared_ptr<Array<Matrix<double>>> ETmats;

  public:
    EmbTrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
                       bool parseflags = false)
        : L2HighOrderFESpace (ama, flags, parseflags)
    {
      name = "EmbTrefftzFESpace";
      type = "embt";
      needs_transform_vec = true;
    }

    shared_ptr<BaseVector>
    SetOp (shared_ptr<SumOfIntegrals> bf, shared_ptr<SumOfIntegrals> lf,
           double eps, shared_ptr<FESpace> test_fes, int tndof);

    void VTransformMR (ElementId ei, const SliceMatrix<double> mat,
                       TRANSFORM_TYPE type) const
    {
      static Timer timer ("EmbTrefftz: MTransform");
      RegionTimer reg (timer);

      size_t nz = ((*ETmats)[ei.Nr ()]).Width ();
      Matrix<double> bla (mat.Height (), mat.Width ());

      if (type == TRANSFORM_MAT_LEFT)
        {
          bla.Rows (0, nz) = Trans ((*ETmats)[ei.Nr ()]) * mat;
          mat = bla;
        }
      if (type == TRANSFORM_MAT_RIGHT)
        {
          bla.Cols (0, nz) = mat * ((*ETmats)[ei.Nr ()]);
          mat = bla;
        }
      if (type == TRANSFORM_MAT_LEFT_RIGHT)
        {
          // mat = mat * Trans(*ETmats[ei.Nr()]);
          // mat = (*ETmats[ei.Nr()]) * mat;
          bla.Cols (0, nz) = mat * ((*ETmats)[ei.Nr ()]);
          mat.Rows (0, nz) = Trans ((*ETmats)[ei.Nr ()]) * bla;
        }
    }

    virtual void VTransformVR (ElementId ei, const SliceVector<double> vec,
                               TRANSFORM_TYPE type) const
    {
      static Timer timer ("EmbTrefftz: VTransform");
      RegionTimer reg (timer);
      size_t nz = ((*ETmats)[ei.Nr ()]).Width ();
      if (type == TRANSFORM_RHS)
        {
          Vector<double> new_vec (vec.Size ());
          new_vec = Trans ((*ETmats)[ei.Nr ()]) * vec;
          vec = new_vec;
        }
      else if (type == TRANSFORM_SOL)
        {
          Vector<double> new_vec (vec.Size ());
          new_vec = ((*ETmats)[ei.Nr ()]) * vec;
          vec = new_vec;
        }
      else
        cout << "transform " << type << " nothing here" << endl;
    }
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEmbTrefftz (py::module m);
#endif // NGS_PYTHON

#endif
