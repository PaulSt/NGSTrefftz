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
  std::tuple<vector<shared_ptr<Matrix<SCAL>>>, shared_ptr<BaseVector>>
  EmbTrefftz (shared_ptr<SumOfIntegrals> bf, shared_ptr<FESpace> fes,
              shared_ptr<SumOfIntegrals> lf, double eps,
              shared_ptr<FESpace> test_fes, int tndof, bool getrange,
              std::map<std::string, Vector<SCAL>> *stats = nullptr);

  template <typename T, typename shrdT>
  class EmbTrefftzFESpace
      : public T //, public std::enable_shared_from_this<EmbTrefftzFESpace>
  {
    vector<shared_ptr<Matrix<double>>> ETmats;
    shrdT fes;
    Array<DofId> all2comp;

  public:
    EmbTrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
                       bool parseflags = false)
        : T (ama, flags, parseflags)
    {
      throw Exception ("Please provide a base fes for the embedding");
    }

    EmbTrefftzFESpace (shrdT afes)
        : T (afes->GetMeshAccess (), afes->GetFlags (), false), fes (afes)
    {
      this->name = "EmbTrefftzFESpace";
      this->type = "embt";
      this->needs_transform_vec = true;
      // this->Update();
      // this->UpdateDofTables();
      // this->UpdateCouplingDofArray();
      // this->FinalizeUpdate();
    }

    shared_ptr<BaseVector>
    SetOp (shared_ptr<SumOfIntegrals> bf, shared_ptr<SumOfIntegrals> lf,
           double eps, shared_ptr<FESpace> test_fes, int tndof);

    void GetDofNrs (ElementId ei, Array<int> &dnums) const override;

    virtual void VTransformMR (ElementId ei, const SliceMatrix<double> mat,
                               TRANSFORM_TYPE type) const override;

    virtual void VTransformVR (ElementId ei, const SliceVector<double> vec,
                               TRANSFORM_TYPE type) const override;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEmbTrefftz (py::module m);
#endif // NGS_PYTHON

#endif
