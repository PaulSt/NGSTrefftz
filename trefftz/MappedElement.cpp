#define FILE_MAPPEDELEMENT_CPP

#include "MappedElement.hpp"

namespace ngfem
{

  void MappedElement ::CalcShape (const BaseMappedIntegrationPoint &mip,
                                  BareSliceVector<Complex> shape) const
  {
    CalcShape (mip,
               SliceVector<double> (ndof, 2 * shape.Dist (),
                                    reinterpret_cast<double *> (&shape (0))));
    SliceVector<double> imag_part (
        ndof, 2 * shape.Dist (), reinterpret_cast<double *> (&shape (0)) + 1);
    imag_part = 0.0;
  }

  void MappedElement ::CalcShape (const BaseMappedIntegrationRule &mir,
                                  SliceMatrix<> shape) const
  {
    for (int i = 0; i < mir.Size (); i++)
      CalcShape (mir[i], shape.Col (i));
  }

  void MappedElement ::CalcShape (const SIMD_BaseMappedIntegrationRule &mir,
                                  BareSliceMatrix<SIMD<double>> shape) const
  {
    throw ExceptionNOSIMD ("SIMD - CalcShape not overloaded");
  }

  double MappedElement ::Evaluate (const BaseMappedIntegrationPoint &mip,
                                   BareSliceVector<double> x) const
  {
    VectorMem<20, double> shape (ndof);
    CalcShape (mip, shape);
    // cout << InnerProduct (shape, x) << endl;
    return InnerProduct (shape, x);
  }

  void MappedElement ::Evaluate (const BaseMappedIntegrationRule &mir,
                                 BareSliceVector<double> coefs,
                                 FlatVector<double> vals) const
  {
    for (size_t i = 0; i < mir.Size (); i++) //.GetNIP(); i++)
      vals (i) = Evaluate (mir[i], coefs);
  }

  void MappedElement ::Evaluate (const SIMD_BaseMappedIntegrationRule &mir,
                                 BareSliceVector<> coefs,
                                 BareVector<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (
        string ("Evaluate (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void MappedElement ::Evaluate (const SIMD_BaseMappedIntegrationRule &mir,
                                 SliceMatrix<> coefs,
                                 BareSliceMatrix<SIMD<double>> values) const
  {
    for (size_t i = 0; i < coefs.Width (); i++)
      Evaluate (mir, coefs.Col (i), values.Row (i));
  }

  void
  MappedElement ::Evaluate (const BaseMappedIntegrationRule &mir,
                            SliceMatrix<> coefs, SliceMatrix<> values) const
  {
    VectorMem<100> shapes (coefs.Height ());
    for (size_t i = 0; i < mir.Size (); i++)
      {
        CalcShape (mir[i], shapes);
        values.Row (i) = Trans (coefs) * shapes;
      }
  }

  void MappedElement ::EvaluateTrans (const BaseMappedIntegrationRule &mir,
                                      FlatVector<double> vals,
                                      BareSliceVector<double> coefs) const
  {
    VectorMem<20, double> shape (ndof);
    coefs.AddSize (ndof) = 0.0;
    for (int i = 0; i < mir.Size (); i++) // GetNIP()
      {
        CalcShape (mir[i], shape);
        coefs.AddSize (ndof) += vals (i) * shape;
      }
  }

  void MappedElement ::AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
                                 BareVector<SIMD<double>> values,
                                 BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD (
        string ("AddTrans (simd) not implemented for class ")
        + typeid (*this).name ());
  }

  void MappedElement ::AddTrans (const SIMD_BaseMappedIntegrationRule &mir,
                                 BareSliceMatrix<SIMD<double>> values,
                                 SliceMatrix<> coefs) const
  {
    for (int i = 0; i < coefs.Width (); i++)
      AddTrans (mir, values.Row (i), coefs.Col (i));
  }

  /*
    void MappedElement ::
    CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir,
                      BareSliceMatrix<SIMD<double>> dshapes) const
    {
      throw ExceptionNOSIMD("SIMD - CalcDShape not overloaded");
    }

    void MappedElement ::
    EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<>
    coefs, BareSliceMatrix<SIMD<double>> values) const
    {
      throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for
    class ")+typeid(*this).name());
    }

    void MappedElement ::
    EvaluateGrad (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs,
    BareSliceMatrix<SIMD<double>> values) const
    {
      throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for
    class ")+typeid(*this).name());
    }

    void MappedElement ::
    AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir,
    BareSliceMatrix<SIMD<double>> values, BareSliceVector<> coefs) const
    {
      throw ExceptionNOSIMD (string("AddGradTrans (simd) not implemented for
    class ")+typeid(*this).name());
    }
  */
}
