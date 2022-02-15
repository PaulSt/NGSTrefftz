#ifndef FILE_PLANEWAVEELEMENT_HPP
#define FILE_PLANEWAVEELEMENT_HPP

#include <fem.hpp>
#include "scalarmappedfe.hpp"

namespace ngfem
{




    template <int D>
    class PlaneWaveElement : public ScalarMappedElement<D>
    {
        private:
            Vec<D> GetDirection(int i) const;
              bool iscomplex = true;
              int conj;

        public:
            PlaneWaveElement(int andof, int aord, ELEMENT_TYPE aeltype, Vec<D> aelcenter = 0, double aelsize = 1, float ac = 1.0, int aconj=1) :
            ScalarMappedElement<D>(andof, aord, Matrix<>(), aeltype, aelcenter, aelsize, ac),
            conj(aconj)
            {;}

            using ScalarMappedElement<D>::CalcMappedDShape;
            using ScalarMappedElement<D>::Evaluate;
            using ScalarMappedElement<D>::EvaluateGrad;
            using ScalarMappedElement<D>::AddGradTrans;

            bool ComplexShapes() const override { return true; }

            Vec<D> EvaluateGrad (const BaseMappedIntegrationPoint & ip, BareSliceVector<> x) const override
            { throw Exception("EvaluateGrad only complex for PW "); }
            void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const override
            { throw Exception("CalcShape only complex for PW "); }
            void CalcShape (const BaseMappedIntegrationRule & mir, SliceMatrix<> shape) const override
            { throw Exception("CalcShape only complex for PW "); }
            void CalcDShape (const BaseMappedIntegrationPoint & mip, BareSliceMatrix<> dshape) const override
            { throw Exception("CalcDShape only complex for PW "); }
            void CalcDShape (const BaseMappedIntegrationRule & mir, BareSliceMatrix<> dshapes) const override
            { throw Exception("CalcDShape only complex for PW "); }
            void CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip, BareSliceMatrix<> hddshape) const override
            { throw Exception("Not implemented for PW "); }
            void CalcShape (const SIMD_BaseMappedIntegrationRule & smir, BareSliceMatrix<SIMD<double>> shape) const override
            {cout << "NO SIMD"<<endl; throw ExceptionNOSIMD("SIMD - CalcShape not overloaded");}
            void CalcDShape (const SIMD_BaseMappedIntegrationRule & smir, BareSliceMatrix<SIMD<double>> dshape) const override
            {cout << "NO SIMD"<<endl; throw ExceptionNOSIMD("SIMD - CalcShape not overloaded");}

            HD NGS_DLL_HEADER virtual void Evaluate (const BaseMappedIntegrationRule & mir, BareSliceVector<Complex> coefs, FlatVector<Complex> vals) const;
            HD NGS_DLL_HEADER virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & mip, BareSliceVector<Complex> x) const;
            void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<Complex> shape) const override;
            void CalcDShape (const BaseMappedIntegrationPoint & mip, BareSliceMatrix<Complex> dshape) const;

            void CalcDShape (const BaseMappedIntegrationRule & mir, BareSliceMatrix<Complex> dshapes) const
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                for (size_t i = 0; i < mir.Size(); i++)
                    CalcDShape (mir[i], dshapes.Cols(i*D,(i+1)*D));
            }
            Vec<D> EvaluateGrad (const BaseMappedIntegrationPoint & ip, BareSliceVector<Complex> x) const
            { throw Exception("OII"); }

            INLINE const FlatMatrixFixWidth<D,Complex> GetDShapeComplex (const BaseMappedIntegrationPoint & mip, LocalHeap & lh) const
            {
                FlatMatrixFixWidth<D,Complex> dshape(this->ndof, lh);
                CalcDShape (mip, dshape);
                return dshape;
            }
            virtual Vec<D,Complex> EvaluateGradComplex (const BaseMappedIntegrationPoint & ip, BareSliceVector<Complex> x) const;

    };



    /// Identity
    template <int D, typename FEL = PlaneWaveElement<D> >
    class DiffOpMappedComplex : public DiffOp<	DiffOpMappedComplex<D, FEL> >
    {
        public:
            enum { DIM = 1 };
            enum { DIM_SPACE = D };
            enum { DIM_ELEMENT = D };
            enum { DIM_DMAT = 1 };
            enum { DIFFORDER = 0 };

            static string Name() { return "Id"; }
            static constexpr bool SUPPORT_PML = true;

            static const FEL & Cast (const FiniteElement & fel)
            { return static_cast<const FEL&> (fel); }

            static void GenerateMatrix (const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
            {
                cout << "Not for complex shapes";
            }
            static void GenerateMatrix (const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    SliceMatrix<Complex,ColMajor> mat, LocalHeap & lh)
            {
                Cast(fel).CalcShape (mip, mat.Row(0));
            }

            template <typename MAT>
            static void GenerateMatrixIR (const FiniteElement & fel,
                    const BaseMappedIntegrationRule & mir,
                    MAT & mat, LocalHeap & lh)
            {
                Cast(fel).CalcShape (mir, Trans(mat));
            }

            template <typename MIP, class TVX, class TVY>
            static void Apply (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                HeapReset hr(lh);
                y = Trans (Cast(fel).GetShape (mip, lh)) * x;
            }

            static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
                    const BareSliceVector<Complex> & x, FlatVector<Complex> & y,
                    LocalHeap & lh)
            {
                y(0) = Cast(fel).EvaluateComplex(mip, x);
            }

            template <class MIR, class TMY>
            static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                    BareSliceVector<double> x, TMY y,
                    LocalHeap & lh)
            {
                cout << "Not for complex shapes";
            }

            template <class MIR>
            static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                    BareSliceVector<Complex> x, SliceMatrix<Complex> y,
                    LocalHeap & lh)
            {
                Cast(fel).Evaluate (mir, x, FlatVector<Complex> (mir.Size(), &y(0,0)));
            }

            template <typename MIP, class TVX, class TVY>
            static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                HeapReset hr(lh);
                y.Range(0,fel.GetNDof()) = Cast(fel).GetShape (mip, lh) * x;
            }

            template <typename MIR, class TVX, class TVY>
            static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                          const TVX & x, TVY & y,
                          LocalHeap & lh)
            {
              y.Range(0,DIM*fel.GetNDof()) = 0.0;
              for (size_t i = 0; i < mir.Size(); i++)
                {
                  HeapReset hr(lh);
                  ApplyTransAdd (fel, mir[i], x.Row(i), y, lh);
                }
            }
            //template <typename MIP, class TVX, class TVY>
            //static void ApplyTransAdd (const FiniteElement & fel, const MIP & mip,
                                       //const TVX & x, TVY & y,
                                       //LocalHeap & lh)
            //{
                //if(fel.ComplexShapes())
                //{
                  //HeapReset hr(lh);
                  //FlatMatrixFixHeight<DIM_DMAT, Complex> mat(DIM*fel.GetNDof(), lh);
                  ////Cast(fel).CalcShape (mip, mat.Row(0));
                  //GenerateMatrix (fel, mip, mat, lh);
                  //y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                //}else
                //{ cout << "Not for complex shapes";}
            //}
    template <typename MIP, class TVX, class TVY>
    static void ApplyTransAdd (const FiniteElement & fel, const MIP & mip,
                               const TVX & x, TVY & y,
                               LocalHeap & lh)
    {
      HeapReset hr(lh);
      typedef typename TVX::TSCAL TSCAL;
      FlatMatrixFixHeight<DIM_DMAT, TSCAL> mat(DIM*fel.GetNDof(), lh);
      //Cast(fel).CalcShape (mip, mat.Row(0));
      GenerateMatrix (fel, mip, mat, lh);
      y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
    }
    };


    /// Gradient operator of dimension D
    template <int D, typename FEL = PlaneWaveElement<D> >
    class DiffOpMappedGradientComplex : public DiffOp<DiffOpMappedGradientComplex<D, FEL> >
    {
        public:
            enum { DIM = 1 };
            enum { DIM_SPACE = D };
            enum { DIM_ELEMENT = D };
            enum { DIM_DMAT = D };
            enum { DIFFORDER = 1 };

            static string Name() { return "grad"; }
            static constexpr bool SUPPORT_PML = true;

            static const FEL & Cast (const FiniteElement & fel)
            {
                return static_cast<const FEL&> (fel);
            }

            template <typename SCALMIP>
            static void GenerateMatrix (const FiniteElement & fel,
                    const MappedIntegrationPoint<D,D,SCALMIP> & mip,
                    SliceMatrix<Complex,ColMajor> mat, LocalHeap & lh)
            {
                Cast(fel).CalcDShape (mip, Trans(mat));
            }

            template <typename SCALMIP, typename MAT>
            static void GenerateMatrix (const FiniteElement & fel,
                    const MappedIntegrationPoint<D,D,SCALMIP> & mip,
                    MAT && mat, LocalHeap & lh)
            {
                HeapReset hr(lh);
                mat = Trans (Cast(fel).GetDShape(mip,lh));
            }

            static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
                    const BareSliceVector<Complex> & x, FlatVector<Complex> && y,
                    LocalHeap & lh)
            {
                Vec<D> hv = Cast(fel).EvaluateGrad(mip, x);
                //y = Trans (mip.GetJacobianInverse()) * hv;
                y = hv;
            }
            template <typename MIP, class TVY>
            static void Apply (const FiniteElement & fel, const MIP & mip,
                    const BareSliceVector<Complex> & x, TVY && y,
                    LocalHeap & lh)
            {
                Vec<D,Complex> hv = Cast(fel).EvaluateGradComplex(mip, x);
                y = hv;
                HeapReset hr(lh);
            }


            template <typename MIP, class TVX, class TVY>
            static void Apply (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY && y,
                    LocalHeap & lh)
            {
                HeapReset hr(lh);
                typedef typename TVX::TSCAL TSCAL;
                Vec<D,TSCAL> hv = Trans (Cast(fel).GetDShape(mip, lh)) * x;
                y = hv;
                //y = Trans (mip.GetJacobianInverse()) * hv;
            }


            //using DiffOp<DiffOpMappedGradient<D, FEL> >::ApplyTans;
            template <typename MIP, class TVX, class TVY>
            static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                typedef typename TVX::TSCAL TSCAL;
                Vec<D,TSCAL> vx = x;
                //auto hv = mip.GetJacobianInverse() * vx;
                y.Range(0,fel.GetNDof()) = Cast(fel).GetDShape(mip,lh) * vx; //* hv;
            }

            template <typename MIP, class TVX, class TVY>
            static void ApplyTransAdd (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                typedef typename MIP::TSCAL TSCAL;

                HeapReset hr(lh);

                FlatMatrixFixHeight<DIM_DMAT, TSCAL> mat(DIM*fel.GetNDof(), lh);
                GenerateMatrix (fel, mip, mat, lh);
                y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
            }
            // for complex shapes in lfi, symbolicintegrator.cpp T_CalcFacetVector -> diffop_impl.hpp ApplyTrans,ApplyTransIR ->
            template <class MIP>
            static void ApplyTransAdd (const FiniteElement & fel,
                    const MIP & mip,
                    FlatVector<Complex> x, BareSliceVector<Complex> y,
                    LocalHeap & lh)
            {
                FlatMatrixFixHeight<DIM_DMAT, Complex> mat(DIM*fel.GetNDof(), lh);
                y.Range(DIM*fel.GetNDof()) += (Cast(fel).GetDShapeComplex(mip,lh)) * x;
                //GenerateMatrix (fel, mip, mat, lh);
                //y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
            }

            template <typename MIR, class TVX, class TVY>
            static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                y.Range(0,DIM*fel.GetNDof()) = 0.0;
                for (size_t i = 0; i < mir.Size(); i++)
                {
                    HeapReset hr(lh);
                    ApplyTransAdd (fel, mir[i], x.Row(i), y, lh);
                }
            }

    };
}
#endif
