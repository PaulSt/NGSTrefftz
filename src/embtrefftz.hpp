#ifndef FILE_SVDTREFFTZ_HPP
#define FILE_SVDTREFFTZ_HPP
#include <comp.hpp>
#include <python_comp.hpp>

#ifndef FILE_INTEGRATORCFHPP
#include <integratorcf.hpp>
#define FILE_INTEGRATORCFHPP
#endif

namespace ngcomp
{

  class TrefftzEmbedding
  {
    shared_ptr<FESpace> fes = nullptr;
    shared_ptr<FESpace> fes_test = nullptr;
    shared_ptr<FESpace> fes_conformity = nullptr;
    const shared_ptr<SumOfIntegrals> top;
    const shared_ptr<SumOfIntegrals> trhs;
    const shared_ptr<SumOfIntegrals> cop;
    const shared_ptr<SumOfIntegrals> crhs;
    shared_ptr<MeshAccess> ma;
    shared_ptr<BitArray> ignoredofs;

    std::variant<size_t, double> ndof_trefftz;
    // shared_ptr<std::map<std::string, Vector<SCAL>>> stats = nullptr;

    shared_ptr<std::map<std::string, Vector<double>>> stats = nullptr;

    /// elmat = (elmat_conforming | elmat_trefftz)
    vector<optional<Matrix<double>>> etmats;
    /// elmat = (elmat_conforming | elmat_trefftz)
    vector<optional<Matrix<Complex>>> etmatsc;
    /// elmat_trefftz.width == local_ndof_trefftz
    /// elmat_conforming.width == elmat.width - local_ndof_trefftz
    ///                        == local_ndof_conforming
    vector<size_t> local_ndofs_trefftz;

    shared_ptr<const BaseVector> psol;

    /// Executes the embedded Trefftz procedure
    /// and stores the results in this class.
    ///
    /// @tparam SCAL either `double` or `Complex`
    template <typename SCAL> void EmbTrefftz ();

  public:
    TrefftzEmbedding (shared_ptr<SumOfIntegrals> _top,
                      shared_ptr<SumOfIntegrals> _trhs,
                      shared_ptr<SumOfIntegrals> _cop,
                      shared_ptr<SumOfIntegrals> _crhs, size_t _ndof_trefftz,
                      double _eps, shared_ptr<FESpace> _fes = nullptr,
                      shared_ptr<FESpace> _fes_test = nullptr,
                      shared_ptr<FESpace> _fes_conformity = nullptr,
                      shared_ptr<BitArray> _ignoredofs = nullptr);

    shared_ptr<BaseVector>
    Embed (const shared_ptr<const BaseVector> tgfu) const;
    shared_ptr<const BaseVector> GetParticularSolution () const;
    shared_ptr<const BaseMatrix> GetEmbedding () const;
    shared_ptr<FESpace> GetFES () const noexcept { return fes; }
    shared_ptr<FESpace> GetFESconf () const noexcept { return fes_conformity; }
    shared_ptr<const BitArray> GetIgnoredDofs () const { return ignoredofs; }
    const vector<optional<Matrix<double>>> &GetEtmats () const noexcept
    {
      return etmats;
    }
    const vector<optional<Matrix<Complex>>> &GetEtmatsC () const noexcept
    {
      return etmatsc;
    }
    const vector<size_t> &GetLocalNodfsTrefftz () const noexcept
    {
      return local_ndofs_trefftz;
    }
  };

  /// Represents the FESpace, that is generated by the embedded Trefftz method.
  /// Use the \ref EmbTrefftzFESpace(shared_ptr<T>) constructor to build the
  /// space from a given (trial) FESpcae, use \ref EmbTrefftzFESpace::SetOp to
  /// specify the operation that the Trefftz space should be build from.
  ///
  /// @tparam T must be a FESpace
  template <typename T>
  class EmbTrefftzFESpace
      : public T //, public std::enable_shared_from_this<EmbTrefftzFESpace>
  {
    shared_ptr<TrefftzEmbedding> emb;
    static_assert (std::is_base_of_v<FESpace, T>, "T must be a FESpace");
    vector<optional<Matrix<double>>> etmats;
    vector<optional<Matrix<Complex>>> etmatsc;
    shared_ptr<T> fes;
    shared_ptr<const FESpace> fes_conformity;
    shared_ptr<const BitArray> ignoredofs;

    /// contains the mapping of Element Number to the associated Dofs
    Table<DofId> elnr_to_dofs;

  public:
    EmbTrefftzFESpace (shared_ptr<MeshAccess> ama, const Flags &flags,
                       bool parseflags = false)
        : T (ama, flags, parseflags)
    {
      throw Exception ("Please construct via an embedding");
    }

    EmbTrefftzFESpace (shared_ptr<TrefftzEmbedding> aemb)
        : T (aemb->GetFES ()->GetMeshAccess (), aemb->GetFES ()->GetFlags (),
             false),
          emb (aemb)
    {
      fes = dynamic_pointer_cast<T> (aemb->GetFES ());
      assert (fes && "fes may not be nullptr");
      this->name = "EmbTrefftzFESpace(" + fes->GetClassName () + ")";
      this->type = "embt";
      this->needs_transform_vec = true;
      this->iscomplex = fes->IsComplex ();
      if constexpr (std::is_same_v<CompoundFESpace, T>)
        for (auto space : fes->Spaces ())
          this->AddSpace (space);

      etmats = emb->GetEtmats ();
      etmatsc = emb->GetEtmatsC ();
      fes_conformity = emb->GetFESconf ();
      ignoredofs = emb->GetIgnoredDofs ();

      adjustDofsAfterSetOp ();
      // this->Update();
      // this->UpdateDofTables();
      // this->UpdateCouplingDofArray();
      // this->FinalizeUpdate();
    }

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;

    virtual void VTransformMR (ElementId ei, const SliceMatrix<double> mat,
                               TRANSFORM_TYPE type) const override;

    virtual void VTransformMC (ElementId ei, const SliceMatrix<Complex> mat,
                               TRANSFORM_TYPE type) const override;

    virtual void VTransformVR (ElementId ei, const SliceVector<double> vec,
                               TRANSFORM_TYPE type) const override;

    virtual void VTransformVC (ElementId ei, const SliceVector<Complex> vec,
                               TRANSFORM_TYPE type) const override;

    virtual string GetClassName () const override;

  private:
    /// adjusts the dofs of the space. Will be called by SetOp.
    void adjustDofsAfterSetOp ();
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEmbTrefftz (py::module m);
#endif // NGS_PYTHON

#endif
