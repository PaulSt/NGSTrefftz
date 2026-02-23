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
  void copyBitArray (const shared_ptr<BitArray> target,
                     const shared_ptr<const BitArray> source);

  template <typename SCAL>
  size_t calcNdofTrefftz (const size_t ndof, const size_t ndof_test,
                          const size_t ndof_conforming,
                          const std::variant<size_t, double> ndof_trefftz,
                          const bool trefftz_op_is_null,
                          const SliceVector<SCAL> singular_values);

  class TrefftzEmbedding
  {
    shared_ptr<FESpace> fes = nullptr;
    shared_ptr<FESpace> fes_test = nullptr;
    shared_ptr<FESpace> fes_conformity = nullptr;
    const shared_ptr<SumOfIntegrals> top;
    const shared_ptr<SumOfIntegrals> trhs;
    const shared_ptr<SumOfIntegrals> cop;
    const shared_ptr<SumOfIntegrals> crhs;
    const shared_ptr<SumOfIntegrals> fes_ip;
    shared_ptr<MeshAccess> ma;
    shared_ptr<BitArray> ignoredofs;

    /// contains the mapping of Element Number to the associated dofs
    /// of the conforming Trefftz space
    Table<DofId> tdof_nrs;

    bool compute_elmat_T_inv = true;

    std::variant<size_t, double> ndof_trefftz;

    shared_ptr<std::map<std::string, Vector<double>>> stats = nullptr;

    /// elmat = (elmat_conforming | elmat_trefftz)
    Array<optional<Matrix<double>>> etmats;
    Array<optional<Matrix<Complex>>> etmatsc;

    /// (pseudo-)inverse matrices of the Trefftz part of etmats
    Array<optional<Matrix<double>>> etmats_trefftz_inv;
    /// (pseudo-)inverse matrices of the Trefftz part of etmatsc
    Array<optional<Matrix<Complex>>> etmatsc_trefftz_inv;

    /// elmat_trefftz.width == local_ndof_trefftz
    /// elmat_conforming.width == elmat.width - local_ndof_trefftz
    ///                        == local_ndof_conforming
    Array<size_t> local_ndofs_trefftz;

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
                      shared_ptr<SumOfIntegrals> _fes_ip = nullptr,
                      shared_ptr<BitArray> _ignoredofs = nullptr,
                      shared_ptr<std::map<std::string, Vector<double>>> _stats
                      = nullptr);

    shared_ptr<BaseVector>
    Embed (const shared_ptr<const BaseVector> tgfu) const;
    shared_ptr<GridFunction>
    Embed (const shared_ptr<const GridFunction> tgfu) const;
    shared_ptr<const BaseVector> GetParticularSolution () const;
    shared_ptr<const BaseVector>
    GetParticularSolution (shared_ptr<SumOfIntegrals> _trhs) const;
    shared_ptr<const BaseVector>
    GetParticularSolution (shared_ptr<const BaseVector> _trhsvec) const;
    shared_ptr<const BaseMatrix> GetEmbedding () const;
    shared_ptr<FESpace> GetFES () const noexcept { return fes; }
    shared_ptr<FESpace> GetFEStest () const noexcept { return fes_test; }
    shared_ptr<FESpace> GetFESconf () const noexcept { return fes_conformity; }
    shared_ptr<const BitArray> GetIgnoredDofs () const { return ignoredofs; }
    const Array<optional<Matrix<double>>> &GetEtmats () const noexcept
    {
      return etmats;
    }
    const Array<optional<Matrix<Complex>>> &GetEtmatsC () const noexcept
    {
      return etmatsc;
    }
    const optional<Matrix<double>> &GetEtmat (size_t i) const
    {
      return etmats[i];
    }
    const optional<Matrix<Complex>> &GetEtmatC (size_t i) const
    {
      return etmatsc[i];
    }
    const Array<size_t> &GetLocalNodfsTrefftz () const noexcept
    {
      return local_ndofs_trefftz;
    }
    const FlatArray<DofId> GetTDofNrs (const size_t elnr) const
    {
      return tdof_nrs[elnr];
    }
  };

  class EmbeddedTrefftzFES : public FESpace
  {
    const shared_ptr<TrefftzEmbedding> emb;
    shared_ptr<FESpace> basefes;

    const FlatArray<optional<Matrix<double>>> etmats;
    const FlatArray<optional<Matrix<Complex>>> etmatsc;
    shared_ptr<FESpace> fes_conformity;
    shared_ptr<const BitArray> ignoredofs;

    /// (pseudo-)inverse matrices of etmats
    mutable Array<optional<Matrix<double>>> etmats_inv;
    mutable Array<optional<Matrix<Complex>>> etmatsc_inv;
    mutable once_flag etmats_inv_computed;

  public:
    EmbeddedTrefftzFES (shared_ptr<MeshAccess> ama, const Flags &flags,
                        bool parseflags = false)
        : FESpace (ama, flags, parseflags)
    {
      throw Exception ("Please construct via a TrefftzEmbedding");
    }

    EmbeddedTrefftzFES (shared_ptr<TrefftzEmbedding> aemb);

    void Update () override;
    void UpdateDofTables () override;
    void UpdateFreeDofs () override;
    void UpdateCouplingDofArray () override;

    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;

    void VTransformMR (ElementId ei, SliceMatrix<double> mat,
                       TRANSFORM_TYPE type) const override;
    void VTransformMC (ElementId ei, SliceMatrix<Complex> mat,
                       TRANSFORM_TYPE type) const override;
    void VTransformVR (ElementId ei, SliceVector<double> vec,
                       TRANSFORM_TYPE type) const override;
    void VTransformVC (ElementId ei, SliceVector<Complex> vec,
                       TRANSFORM_TYPE type) const override;

    string GetClassName () const override;

    shared_ptr<TrefftzEmbedding> GetEmbedding () const noexcept { return emb; }
    shared_ptr<FESpace> GetBaseFESpace () const noexcept { return basefes; }

    virtual ProxyNode MakeProxyFunction (
        bool testfunction,
        const function<shared_ptr<ProxyFunction> (shared_ptr<ProxyFunction>)>
            &addblock) const override;

    virtual SymbolTable<shared_ptr<DifferentialOperator>>
    GetAdditionalEvaluators () const override
    {
      return basefes->GetAdditionalEvaluators ();
    }

  private:
    optional<FlatMatrix<double>> GetEtmatInv (size_t idx) const;
    optional<FlatMatrix<Complex>> GetEtmatCInv (size_t idx) const;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEmbTrefftz (py::module m);
#endif // NGS_PYTHON

#endif
