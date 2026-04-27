#ifndef FILE_PLANEWAVEELEMENT_HPP
#define FILE_PLANEWAVEELEMENT_HPP

#include <fem.hpp>
#include "ngsttd.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{
  template <int D> class PlaneWaveElement : public BaseScalarMappedElement
  {
  private:
    Vec<D> GetDirection (int i) const;

    ELEMENT_TYPE eltype;
    Vec<D> shift;

    double elsize;
    double c;
    int conj;

  public:
    PlaneWaveElement (int andof, int aord, ELEMENT_TYPE aeltype,
                      Vec<D> ashift = 0, double aelsize = 1, double ac = 1.0,
                      int aconj = 1)
        : BaseScalarMappedElement (andof, aord), eltype (aeltype),
          shift (ashift), elsize (aelsize), c (ac), conj (aconj)
    {
      ;
    }

    bool ComplexShapes () const override { return true; }

    virtual ELEMENT_TYPE ElementType () const override { return eltype; }

    using BaseScalarMappedElement::AddGradTrans;
    using BaseScalarMappedElement::CalcDShape;
    using BaseScalarMappedElement::CalcShape;
    using BaseScalarMappedElement::Evaluate;
    using BaseScalarMappedElement::EvaluateGrad;

    Vec<D>
    EvaluateGrad (const BaseMappedIntegrationPoint &, BareSliceVector<>) const
    {
      throw Exception ("EvaluateGrad only complex for PW ");
    }
    void CalcShape (const BaseMappedIntegrationPoint &,
                    BareSliceVector<>) const override
    {
      throw Exception ("CalcShape only complex for PW ");
    }
    void CalcShape (const BaseMappedIntegrationRule &,
                    BareSliceMatrix<>) const override
    {
      throw Exception ("CalcShape only complex for PW ");
    }
    void CalcDShape (const BaseMappedIntegrationPoint &,
                     BareSliceMatrix<>) const override
    {
      throw Exception ("CalcDShape only complex for PW ");
    }

    void
    CalcDShape (const BaseMappedIntegrationRule &, BareSliceMatrix<>) const
    {
      throw Exception ("CalcDShape only complex for PW ");
    }
    void CalcMappedDDShape (const BaseMappedIntegrationPoint &,
                            BareSliceMatrix<>) const
    {
      throw Exception ("Not implemented for PW ");
    }
    void CalcShape (const SIMD_BaseMappedIntegrationRule &,
                    BareSliceMatrix<SIMD<double>>) const override
    {
      throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
    }
    void CalcDShape (const SIMD_BaseMappedIntegrationRule &,
                     BareSliceMatrix<SIMD<double>>) const
    {
      throw ExceptionNOSIMD ("SIMD - CalcDShape not overloaded");
    }

    NGST_DLL virtual void
    Evaluate (const BaseMappedIntegrationRule &mir,
              BareSliceVector<Complex> coefs, FlatVector<Complex> vals) const;
    NGST_DLL virtual Complex
    EvaluateComplex (const BaseMappedIntegrationPoint &mip,
                     BareSliceVector<Complex> x) const;
    void CalcShape (const BaseMappedIntegrationPoint &mip,
                    BareSliceVector<Complex> shape) const override;
    void CalcDShape (const BaseMappedIntegrationPoint &mip,
                     BareSliceMatrix<Complex> dshape) const;

    void CalcDShape (const BaseMappedIntegrationRule &mir,
                     BareSliceMatrix<Complex> dshapes) const
    {
      for (size_t i = 0; i < mir.Size (); i++)
        CalcDShape (mir[i], dshapes.Cols (i * D, (i + 1) * D));
    }
    Vec<D> EvaluateGrad (const BaseMappedIntegrationPoint &,
                         BareSliceVector<Complex>) const
    {
      throw Exception ("OII");
    }

    INLINE const FlatMatrixFixWidth<D, Complex>
    GetDShapeComplex (const BaseMappedIntegrationPoint &mip,
                      LocalHeap &lh) const
    {
      FlatMatrixFixWidth<D, Complex> dshape (this->ndof, lh);
      CalcDShape (mip, dshape);
      return dshape;
    }

    virtual Vec<D, Complex>
    EvaluateGradComplex (const BaseMappedIntegrationPoint &ip,
                         BareSliceVector<Complex> x) const;
  };

  // --------------------------------------------------------------------------
  // Direct complex-valued value operator for PlaneWaveElement
  // --------------------------------------------------------------------------
  template <int D> class PlaneWaveValueOperator : public DifferentialOperator
  {
  public:
    PlaneWaveValueOperator () : DifferentialOperator (1, 1, VOL, 0) { ; }

    string Name () const override { return "pw-id"; }

  private:
    static const PlaneWaveElement<D> &Cast (const FiniteElement &fel)
    {
      return static_cast<const PlaneWaveElement<D> &> (fel);
    }

  public:
    void CalcMatrix (const FiniteElement &fel,
                     const BaseMappedIntegrationPoint &mip,
                     BareSliceMatrix<Complex, ColMajor> mat,
                     LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      FlatVector<Complex> shape (ndof, lh);
      Cast (fel).CalcShape (mip, shape);

      for (size_t j = 0; j < ndof; ++j)
        mat (0, j) = shape (j);
    }

    void
    CalcMatrix (const FiniteElement &fel, const BaseMappedIntegrationRule &mir,
                BareSliceMatrix<Complex, ColMajor> mat,
                LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      for (size_t ip = 0; ip < mir.Size (); ++ip)
        {
          HeapReset hr (lh);

          FlatVector<Complex> shape (ndof, lh);
          Cast (fel).CalcShape (mir[ip], shape);

          for (size_t j = 0; j < ndof; ++j)
            mat (ip, j) = shape (j);
        }
    }

    void
    Apply (const FiniteElement &fel, const BaseMappedIntegrationPoint &mip,
           BareSliceVector<Complex> x, FlatVector<Complex> value,
           LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      FlatVector<Complex> shape (ndof, lh);
      Cast (fel).CalcShape (mip, shape);

      Complex sum = 0.0;
      for (size_t j = 0; j < ndof; ++j)
        sum += shape (j) * x (j);

      value (0) = sum;
    }

    void Apply (const FiniteElement &fel, const BaseMappedIntegrationRule &mir,
                BareSliceVector<Complex> x, BareSliceMatrix<Complex> values,
                LocalHeap &lh) const override
    {
      auto vals = values.AddSize (mir.Size (), 1);

      for (size_t ip = 0; ip < mir.Size (); ++ip)
        {
          HeapReset hr (lh);

          FlatVector<Complex> shape (fel.GetNDof (), lh);
          Cast (fel).CalcShape (mir[ip], shape);

          Complex sum = 0.0;
          for (int j = 0; j < fel.GetNDof (); ++j)
            sum += shape (j) * x (j);

          vals (ip, 0) = sum;
        }
    }

    void ApplyTrans (const FiniteElement &fel,
                     const BaseMappedIntegrationPoint &mip,
                     FlatVector<Complex> value, BareSliceVector<Complex> x,
                     LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      auto xr = x.Range (0, ndof);
      xr = 0.0;

      FlatVector<Complex> shape (ndof, lh);
      Cast (fel).CalcShape (mip, shape);

      // This is transpose, not conjugate transpose.
      for (size_t j = 0; j < ndof; ++j)
        xr (j) += value (0) * shape (j);
    }

    void
    ApplyTrans (const FiniteElement &fel, const BaseMappedIntegrationRule &mir,
                FlatMatrix<Complex> values, BareSliceVector<Complex> x,
                LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      auto xr = x.Range (0, ndof);
      xr = 0.0;

      for (size_t ip = 0; ip < mir.Size (); ++ip)
        {
          HeapReset hr (lh);

          FlatVector<Complex> shape (ndof, lh);
          Cast (fel).CalcShape (mir[ip], shape);

          for (size_t j = 0; j < ndof; ++j)
            xr (j) += values (ip, 0) * shape (j);
        }
    }
  };

  // --------------------------------------------------------------------------
  // Direct complex-valued gradient operator for PlaneWaveElement
  // --------------------------------------------------------------------------
  template <int D>
  class PlaneWaveGradientOperator : public DifferentialOperator
  {
  private:
    static const PlaneWaveElement<D> &Cast (const FiniteElement &fel)
    {
      return static_cast<const PlaneWaveElement<D> &> (fel);
    }

  public:
    PlaneWaveGradientOperator () : DifferentialOperator (D, 1, VOL, 1) { ; }

    string Name () const override { return "pw-grad"; }

    void CalcMatrix (const FiniteElement &fel,
                     const BaseMappedIntegrationPoint &mip,
                     BareSliceMatrix<Complex, ColMajor> mat,
                     LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      FlatMatrixFixWidth<D, Complex> dshape (ndof, lh);
      Cast (fel).CalcDShape (mip, dshape);

      // dshape(j,d) = d phi_j / d x_d.
      // B-matrix convention: rows are components, columns are local dofs.
      for (int d = 0; d < D; ++d)
        for (size_t j = 0; j < ndof; ++j)
          mat (d, j) = dshape (j, d);
    }

    void
    CalcMatrix (const FiniteElement &fel, const BaseMappedIntegrationRule &mir,
                BareSliceMatrix<Complex, ColMajor> mat,
                LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      for (size_t ip = 0; ip < mir.Size (); ++ip)
        {
          HeapReset hr (lh);

          FlatMatrixFixWidth<D, Complex> dshape (ndof, lh);
          Cast (fel).CalcDShape (mir[ip], dshape);

          for (int d = 0; d < D; ++d)
            for (size_t j = 0; j < ndof; ++j)
              mat (ip * D + d, j) = dshape (j, d);
        }
    }

    void
    Apply (const FiniteElement &fel, const BaseMappedIntegrationPoint &mip,
           BareSliceVector<Complex> x, FlatVector<Complex> grad,
           LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      FlatMatrixFixWidth<D, Complex> dshape (ndof, lh);
      Cast (fel).CalcDShape (mip, dshape);

      for (int d = 0; d < D; ++d)
        {
          Complex sum = 0.0;
          for (size_t j = 0; j < ndof; ++j)
            sum += dshape (j, d) * x (j);

          grad (d) = sum;
        }
    }

    void Apply (const FiniteElement &fel, const BaseMappedIntegrationRule &mir,
                BareSliceVector<Complex> x, BareSliceMatrix<Complex> grads,
                LocalHeap &lh) const override
    {
      auto gradmat = grads.AddSize (mir.Size (), D);

      for (size_t ip = 0; ip < mir.Size (); ++ip)
        {
          HeapReset hr (lh);

          FlatMatrixFixWidth<D, Complex> dshape (fel.GetNDof (), lh);
          Cast (fel).CalcDShape (mir[ip], dshape);

          for (int d = 0; d < D; ++d)
            {
              Complex sum = 0.0;
              for (int j = 0; j < fel.GetNDof (); ++j)
                sum += dshape (j, d) * x (j);

              gradmat (ip, d) = sum;
            }
        }
    }

    void ApplyTrans (const FiniteElement &fel,
                     const BaseMappedIntegrationPoint &mip,
                     FlatVector<Complex> grad, BareSliceVector<Complex> x,
                     LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      auto xr = x.Range (0, ndof);
      xr = 0.0;

      FlatMatrixFixWidth<D, Complex> dshape (ndof, lh);
      Cast (fel).CalcDShape (mip, dshape);

      // This is transpose, not conjugate transpose.
      for (size_t j = 0; j < ndof; ++j)
        {
          Complex sum = 0.0;
          for (int d = 0; d < D; ++d)
            sum += dshape (j, d) * grad (d);

          xr (j) += sum;
        }
    }

    void
    ApplyTrans (const FiniteElement &fel, const BaseMappedIntegrationRule &mir,
                FlatMatrix<Complex> grads, BareSliceVector<Complex> x,
                LocalHeap &lh) const override
    {
      const size_t ndof = fel.GetNDof ();

      auto xr = x.Range (0, ndof);
      xr = 0.0;

      for (size_t ip = 0; ip < mir.Size (); ++ip)
        {
          HeapReset hr (lh);

          FlatMatrixFixWidth<D, Complex> dshape (ndof, lh);
          Cast (fel).CalcDShape (mir[ip], dshape);

          for (size_t j = 0; j < ndof; ++j)
            {
              Complex sum = 0.0;
              for (int d = 0; d < D; ++d)
                sum += dshape (j, d) * grads (ip, d);

              xr (j) += sum;
            }
        }
    }
  };

}
#endif
