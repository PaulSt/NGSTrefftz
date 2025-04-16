#include "planewavefe.hpp"

namespace ngfem
{

  template <> Vec<2> PlaneWaveElement<2>::GetDirection (int i) const
  {
    return { cos (2.0 * M_PI * i / this->ndof),
             sin (2.0 * M_PI * i / this->ndof) };
  }

  template <> Vec<3> PlaneWaveElement<3>::GetDirection (int i) const
  {
    int p = this->ndof;
    int q = sqrt (p) - 1;
    int l = 0;
    while (i > 0 && (i / pow (l + 1, 2) >= 1))
      l++;
    int m = i - l * l;

    double theta
        = l % 2 ? M_PI / 2 * l / (q + 1) : M_PI - M_PI / 2 * (l + 1) / (q + 1);
    double phi = l == 0 ? M_PI : M_PI * 2 * m / (2 * l + 1);

    return { sin (theta) * cos (phi), sin (theta) * sin (phi), cos (theta) };
  }

  // template<>
  // Matrix<Complex> PlaneWaveElement<2> :: StableMat() const
  //{
  // Matrix<Complex> M(this->ndof);
  // for (int i=0; i<this->ndof; ++i)
  // for (int j=0; j<this->ndof; ++j)
  //{
  // M(i,j)=exp(-1i*(2.0*M_PI*i/this->ndof)*(j-this->ndof/2.0)*conj);
  //}

  // CalcInverse(M,INVERSE_LIB::INV_LAPACK);
  // return M;
  //}

  template <int D>
  void PlaneWaveElement<D>::CalcShape (const BaseMappedIntegrationPoint &mip,
                                       BareSliceVector<Complex> shape) const
  {
    Vec<D> cpoint = mip.GetPoint ();
    cpoint -= this->shift;

    for (int i = 0; i < this->ndof; ++i)
      {
        Vec<D> dir = GetDirection (i);
        shape (i) = exp (Complex (0, InnerProduct (cpoint, dir) * conj));
      }
  }

  template <int D>
  void PlaneWaveElement<D>::CalcDShape (const BaseMappedIntegrationPoint &mip,
                                        BareSliceMatrix<Complex> dshape) const
  {
    Vec<D> cpoint = mip.GetPoint ();
    cpoint -= this->shift;

    for (int d = 0; d < D; d++)
      for (int i = 0; i < this->ndof; ++i)
        {
          Vec<D> dir = GetDirection (i);
          dshape (i, d)
              = Complex (0, dir[d] * conj)
                * exp (Complex (0, InnerProduct (cpoint, dir) * conj));
        }
  }

  template <int D>
  Complex
  PlaneWaveElement<D>::EvaluateComplex (const BaseMappedIntegrationPoint &mip,
                                        BareSliceVector<Complex> x) const
  {
    VectorMem<20, Complex> shape (this->ndof);
    CalcShape (mip, shape);
    return InnerProduct (shape, x);
  }

  template <int D>
  void PlaneWaveElement<D>::Evaluate (const BaseMappedIntegrationRule &mir,
                                      BareSliceVector<Complex> coefs,
                                      FlatVector<Complex> vals) const
  {
    for (size_t i = 0; i < mir.Size (); i++) //.GetNIP(); i++)
      vals (i) = EvaluateComplex (mir[i], coefs);
  }

  template <int D>
  Vec<D, Complex> PlaneWaveElement<D>::EvaluateGradComplex (
      const BaseMappedIntegrationPoint &ip, BareSliceVector<Complex> x) const
  {
    MatrixFixWidth<D, Complex> dshape (this->ndof);
    CalcDShape (ip, dshape);
    Vec<D, Complex> grad = Trans (dshape) * x;
    return grad;
  }

  template class PlaneWaveElement<2>;
  template class PlaneWaveElement<3>;
}
