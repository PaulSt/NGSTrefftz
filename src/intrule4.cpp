#include <fem.hpp>

namespace ngfem
{
  template <int DIM_ELEMENT, int DIM_SPACE>
  SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> ::
  SIMD_MappedIntegrationRule (const SIMD_IntegrationRule & ir,
                              const ElementTransformation & aeltrans,
                              Allocator & lh)
    : SIMD_BaseMappedIntegrationRule (ir, aeltrans), mips(ir.Size(), lh)
  {
    dim_element = DIM_ELEMENT;
    dim_space = DIM_SPACE;
    baseip = (char*)(void*)(SIMD<BaseMappedIntegrationPoint>*)(&mips[0]);
    incr = sizeof (SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>);

    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    FlatArray<SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>> hmips = mips;
    for (size_t i = 0; i < hir.Size(); i++)
      new (&hmips[i]) SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>> (hir[i], eltrans, -1);

    new (&points) BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                 &mips[0].Point()(0),
                                                 DummySize(mips.Size(), DIM_SPACE));

    new (&normals) BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                  &mips[0].NV()(0),
                                                  DummySize(mips.Size(), DIM_SPACE));

    eltrans.CalcMultiPointJacobian (ir, *this);

    //if (ir.Size())
      //if (ir[0].VB() != VOL)
        //ComputeNormalsAndMeasure (eltrans.GetElementType(), ir[0].FacetNr());
  }

  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> ::
  ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr)
  {throw Exception("Not implemented for SIMD_MappedIntegrationRule<3,4>");}

  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> ::
  TransformGradient (BareSliceMatrix<SIMD<double>> grad) const
  {throw Exception("Not implemented for SIMD_MappedIntegrationRule<3,4>");}

  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> ::
  TransformGradientTrans (BareSliceMatrix<SIMD<double>> grad) const
  {throw Exception("Not implemented for SIMD_MappedIntegrationRule<3,4>");}

  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> :: Print (ostream & ost) const
  {
    ost << "simd-mir, size = " << mips.Size() << endl;
    for (size_t i = 0; i < mips.Size(); i++)
      mips[i].Print(ost);
  }

  template class SIMD_MappedIntegrationRule<3,4>;
}
