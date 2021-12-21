#ifndef FILE_DIFFOPMAPPED
#define FILE_DIFFOPMAPPED

namespace ngfem
{

    /*
       realizations of bdb integrators for many equations.
       The differential operators provide the B-matrix,
       the DMatOps provide the coefficient tensors
       */

    /// Identity
    template <int D, typename FEL = ScalarMappedElement<D> >
    class DiffOpMapped : public DiffOp<	DiffOpMapped<D, FEL> >
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

            //template <typename MIP, typename MAT>
            //static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
                    //MAT && mat, LocalHeap & lh)
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //cout << typeid(mat).name() << endl;
                //HeapReset hr(lh);
                //mat.Row(0) = Cast(fel).GetShape(mip, lh);
            //}

            //static void GenerateMatrix (const FiniteElement & fel,
                    //const BaseMappedIntegrationPoint & mip,
                    //FlatMatrixFixHeight<1> & mat, LocalHeap & lh)
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //Cast(fel).CalcShape (mip, mat.Row(0)); // FlatVector<> (fel.GetNDof(), &mat(0,0)));
            //}
            //static void GenerateMatrix (const FiniteElement & fel,
                    //const BaseMappedIntegrationPoint & mip,
                    //FlatMatrixFixHeight<1,Complex> & mat, LocalHeap & lh)
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //Cast(fel).CalcShape (mip, mat.Row(0)); // FlatVector<> (fel.GetNDof(), &mat(0,0)));
            //}

            static void GenerateMatrix (const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcShape (mip, mat.Row(0));
            }
            static void GenerateMatrix (const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    SliceMatrix<Complex,ColMajor> mat, LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcShape (mip, mat.Row(0));
            }

            // using DiffOp<DiffOpId<D, FEL> > :: GenerateMatrixIR;
            template <typename MAT>
            static void GenerateMatrixIR (const FiniteElement & fel,
                    const BaseMappedIntegrationRule & mir,
                    MAT & mat, LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcShape (mir, Trans(mat));
            }

            static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                    const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> mat)
            {
                Cast(fel).CalcShape (mir, mat);
            }

            template <typename MIP, class TVX, class TVY>
            static void Apply (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
      cout << typeid(x).name() << endl;
      cout << typeid(y).name() << endl;
                HeapReset hr(lh);
                y = Trans (Cast(fel).GetShape (mip, lh)) * x;
            }

            static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
                    const FlatVector<double> & x, FlatVector<double> & y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                y(0) = Cast(fel).Evaluate(mip, x);
            }

            static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
                    const BareSliceVector<Complex> & x, FlatVector<Complex> & y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                y(0) = Cast(fel).EvaluateComplex(mip, x);
            }

            // using DiffOp<DiffOpId<D, FEL> >::ApplyIR;
            template <class MIR, class TMY>
            static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                    BareSliceVector<double> x, TMY y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).Evaluate (mir, x, FlatVector<> (mir.Size(), &y(0,0)));
            }

            template <class MIR>
            static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                    BareSliceVector<Complex> x, SliceMatrix<Complex> y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                //Cast(fel).Evaluate (mir,
                        //SliceMatrix<double> (fel.GetNDof(), 2, 2, reinterpret_cast<double*> (&x(0))),
                        //SliceMatrix<double> (mir.Size(), 2, 2, reinterpret_cast<double*> (&y(0))));
                Cast(fel).Evaluate (mir, x, FlatVector<Complex> (mir.Size(), &y(0,0)));
            }

            // using ApplySIMDIR;
            using DiffOp<DiffOpMapped<D, FEL> >::ApplySIMDIR;
            static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).Evaluate (mir, x, y.Row(0));
            }

            template <typename MIP, class TVX, class TVY>
            static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                HeapReset hr(lh);
                y.Range(0,fel.GetNDof()) = Cast(fel).GetShape (mip, lh) * x;
            }

            // using DiffOp<DiffOpId<D, FEL> >::ApplyTransIR;
            template <class MIR>
            static void ApplyTransIR (const FiniteElement & fel,
                    const MIR & mir,
                    FlatMatrix<double> x, BareSliceVector<double> y,
                    LocalHeap & lh)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).EvaluateTrans (mir, FlatVector<> (mir.Size(), &x(0,0)), y);
            }

            //template <class MIR>
            //static void ApplyTransIR (const FiniteElement & fel,
                    //const MIR & mir,
                    //FlatMatrix<Complex> x, BareSliceVector<Complex> y,
                    //LocalHeap & lh)
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //DiffOp<DiffOpMapped<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);
            //}
    template <typename MIR, class TVX, class TVY>
    static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                  const TVX & x, TVY & y,
                  LocalHeap & lh)
    {
        cout << __FILE__<<" "<<__LINE__<<endl;
      y.Range(0,DIM*fel.GetNDof()) = 0.0;
      for (size_t i = 0; i < mir.Size(); i++)
        {
          HeapReset hr(lh);
          ApplyTransAdd (fel, mir[i], x.Row(i), y, lh);
        }
    }
    template <typename MIP, class TVX, class TVY>
    static void ApplyTransAdd (const FiniteElement & fel, const MIP & mip,
                               const TVX & x, TVY & y,
                               LocalHeap & lh)
    {
        cout << __FILE__<<" "<<__LINE__<<endl;
      //cout << __FILE__<<" "<<__LINE__<<endl;
      //cout << typeid(x).name() << endl;
      //cout << typeid(y).name() << endl;
                if(fel.ComplexShapes())
                {
        cout << __FILE__<<" "<<__LINE__<<endl;
      HeapReset hr(lh);
      FlatMatrixFixHeight<DIM_DMAT, Complex> mat(DIM*fel.GetNDof(), lh);
      //Cast(fel).CalcShape (mip, mat.Row(0));
      GenerateMatrix (fel, mip, mat, lh);
      y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                }
                else
                {
      typedef typename MIP::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DIM_DMAT, TSCAL> mat(DIM*fel.GetNDof(), lh);
      GenerateMatrix (fel, mip, mat, lh);
      y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                }
    }
    //using DiffOp<DiffOpMappedGradient<D, FEL> >::AddTransIR;
    //using DiffOp<DiffOpMappedGradient<D, FEL> >::AddTransAdd;
            // for complex shapes in lfi, symbolicintegrator.cpp T_CalcFacetVector -> diffop_impl.hpp ApplyTrans,ApplyTransIR ->
            //template <class MIP>
            //static void ApplyTransAdd (const FiniteElement & fel,
                    //const MIP & mip,
                    //FlatVector<Complex> x, BareSliceVector<Complex> y,
                    //LocalHeap & lh)
            //{
        //cout << __FILE__<<" "<<__LINE__<<endl;
                //if(fel.ComplexShapes())
                //{
      //FlatMatrixFixHeight<DIM_DMAT, Complex> mat(DIM*fel.GetNDof(), lh);
      ////Cast(fel).CalcShape (mip, mat.Row(0));
      //GenerateMatrix (fel, mip, mat, lh);
      //y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                //}
                //else
                //{
      //FlatMatrixFixHeight<DIM_DMAT, double> mat(DIM*fel.GetNDof(), lh);
      //GenerateMatrix (fel, mip, mat, lh);
      //y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                //}
            //}


            using DiffOp<DiffOpMapped<D, FEL> >::AddTransSIMDIR;
            static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
            {
        cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).AddTrans (mir, y.Row(0), x);
            }

    };



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Identity on boundary
    template <int D, typename FEL = ScalarMappedElement<D-1> >
    class DiffOpMappedBoundary : public DiffOp<DiffOpMappedBoundary<D, FEL> >
    {
        public:
            enum { DIM = 1 };
            enum { DIM_SPACE = D };
            enum { DIM_ELEMENT = D-1 };
            enum { DIM_DMAT = 1 };
            enum { DIFFORDER = 0 };

            static const FEL & Cast (const FiniteElement & fel)
            { return static_cast<const FEL&> (fel); }

            template <typename AFEL, typename MIP, typename MAT>
            static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                    MAT & mat, LocalHeap & lh)
            {
                HeapReset hr(lh);
                mat.Row(0) = static_cast<const FEL&>(fel).GetShape(mip, lh);
            }

            template <typename AFEL, typename MIP, class TVX, class TVY>
            static void Apply (const AFEL & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                HeapReset hr(lh);
                y = Trans (static_cast<const FEL&>(fel).GetShape (mip, lh)) * x;
                // y(0) = InnerProduct (x, static_cast<const FEL&>(fel).GetShape (mip.IP(), lh));
            }

            static void Apply (const ScalarFiniteElement<D-1> & fel, const MappedIntegrationPoint<D-1,D> & mip,
                    const FlatVector<double> & x, FlatVector<double> & y,
                    LocalHeap & lh)
            {
                y(0) = static_cast<const FEL&>(fel).Evaluate(mip, x);
            }

            template <typename AFEL, typename MIP, class TVX, class TVY>
            static void ApplyTrans (const AFEL & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                HeapReset hr(lh);
                y = static_cast<const FEL&>(fel).GetShape (mip, lh) * x;
            }


            // using DiffOp<DiffOpIdBoundary<D, FEL> >::ApplyTransIR;
            template <class MIR>
            static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                    FlatMatrix<double> x, FlatVector<double> y,
                    LocalHeap & lh)
            {
                // static Timer t("applytransir - bnd");
                // RegionTimer reg(t);

                static_cast<const FEL&>(fel).EvaluateTrans (mir, FlatVector<> (mir.Size(), &x(0,0)), y);
            }

            template <class MIR>
            static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                    FlatMatrix<Complex> x, FlatVector<Complex> y,
                    LocalHeap & lh)
            {
                DiffOp<DiffOpMappedBoundary<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);
                // static_cast<const FEL&>(fel).
                // EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
            }

            using DiffOp<DiffOpMappedBoundary<D, FEL> >::ApplySIMDIR;
            static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
            {
                Cast(fel).Evaluate (mir, x, y.Row(0));
            }

            using DiffOp<DiffOpMappedBoundary<D, FEL> >::AddTransSIMDIR;
            static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
            {
                Cast(fel).AddTrans (mir, y.Row(0), x);
            }
    };



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    /// Gradient operator of dimension D
    template <int D, typename FEL = ScalarMappedElement<D> >
    class DiffOpMappedGradient : public DiffOp<DiffOpMappedGradient<D, FEL> >
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

            static void GenerateMatrix (const FiniteElement & fel,
                    const MappedIntegrationPoint<D,D> & mip,
                    SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcMappedDShape (mip, Trans(mat));
            }

            template <typename SCALMIP>
            static void GenerateMatrix (const FiniteElement & fel,
                    const MappedIntegrationPoint<D,D,SCALMIP> & mip,
                    SliceMatrix<Complex,ColMajor> mat, LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcDShape (mip, Trans(mat));
                //cout << mat << endl;
            }

            template <typename SCALMIP, typename MAT>
            static void GenerateMatrix (const FiniteElement & fel,
                    const MappedIntegrationPoint<D,D,SCALMIP> & mip,
                    MAT && mat, LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                cout << typeid(mat).name() << endl;
                cout << typeid(mip).name() << endl;
                ////cout << (mat).IsComplex() << endl;
                //cout << (mip).IsComplex() << endl;
                //cout << (fel.ComplexShapes());
                HeapReset hr(lh);
                //if(fel.ComplexShapes())
                    //mat = Trans (Cast(fel).GetDShapeComplex(mip,lh));
                //else
                    mat = Trans (Cast(fel).GetDShape(mip,lh));
            }

            static void GenerateMatrixIR (const FiniteElement & fel,
                    const MappedIntegrationRule<D,D> & mir,
                    SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcMappedDShape (mir, Trans(mat));
            }

            static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                    const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> mat)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).CalcMappedDShape (mir, mat);
            }

            ///

            //template <class TVY>
            //static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
                    //const FlatVector<> & x, TVY && y,
                    //LocalHeap & lh)
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //Vec<D> hv = Cast(fel).EvaluateGrad(mip, x);
                ////y = Trans (mip.GetJacobianInverse()) * hv;
                //y = hv;
            //}

            static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
                    const BareSliceVector<Complex> & x, FlatVector<Complex> && y,
                    LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Vec<D> hv = Cast(fel).EvaluateGrad(mip, x);
                //y = Trans (mip.GetJacobianInverse()) * hv;
                y = hv;
            }
            template <typename MIP, class TVY>
            static void Apply (const FiniteElement & fel, const MIP & mip,
                    const BareSliceVector<Complex> & x, TVY && y,
                    LocalHeap & lh)
            {
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //cout << typeid(x).name() << endl;
                //cout << typeid(y).name() << endl;
                Vec<D,Complex> hv = Cast(fel).EvaluateGradComplex(mip, x);
                y = hv;
                HeapReset hr(lh);
            }


            template <typename MIP, class TVX, class TVY>
            static void Apply (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY && y,
                    LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                cout << typeid(x).name() << endl;
                cout << typeid(y).name() << endl;
                HeapReset hr(lh);
                typedef typename TVX::TSCAL TSCAL;
                Vec<D,TSCAL> hv = Trans (Cast(fel).GetDShape(mip, lh)) * x;
                y = hv;
                //y = Trans (mip.GetJacobianInverse()) * hv;
            }

            using DiffOp<DiffOpMappedGradient<D, FEL> >::ApplyIR;

            template <class MIR>
            static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                    const FlatVector<double> x, FlatMatrixFixWidth<D,double> y,
                    LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                FlatMatrixFixWidth<D> grad(mir.Size(), &y(0));
                Cast(fel).EvaluateGrad (mir, x, grad);
                /*
                   for (int i = 0; i < mir.Size(); i++)
                   {
                   Vec<D> hv = grad.Row(i);
                   grad.Row(i) = Trans (mir[i].GetJacobianInverse()) * hv;
                   }
                   */
            }

            using DiffOp<DiffOpMappedGradient<D, FEL> >::ApplySIMDIR;
            static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).EvaluateGrad (mir, x, y);
            }


            ///
            //void ApplyTrans (const FiniteElement & fel,
                          //const BaseMappedIntegrationPoint & mip,
                          //FlatVector<Complex> x,
                          //BareSliceVector<Complex> y,
                          //LocalHeap & lh) const override
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //cout << "HELLO";
                ////y.Range(0,fel.GetNDof()) = Cast(fel).GetDShape(mip,lh) * x;
            //}
            //void ApplyTrans (const FiniteElement & fel,
                //const BaseMappedIntegrationRule & mir,
                //FlatMatrix<Complex> flux,
                //BareSliceVector<Complex> x,
                //LocalHeap & lh) const override
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                //cout << "HELLO";
            //}

            //using DiffOp<DiffOpMappedGradient<D, FEL> >::ApplyTans;
            template <typename MIP, class TVX, class TVY>
            static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
                    const TVX & x, TVY & y,
                    LocalHeap & lh)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
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
      cout << __FILE__<<" "<<__LINE__<<endl;
      //cout << typeid(x).name() << endl;
      //cout << typeid(y).name() << endl;
      typedef typename MIP::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DIM_DMAT, TSCAL> mat(DIM*fel.GetNDof(), lh);
      GenerateMatrix (fel, mip, mat, lh);
      y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
    }
    //using DiffOp<DiffOpMappedGradient<D, FEL> >::AddTransIR;
    //using DiffOp<DiffOpMappedGradient<D, FEL> >::AddTransAdd;
            // for complex shapes in lfi, symbolicintegrator.cpp T_CalcFacetVector -> diffop_impl.hpp ApplyTrans,ApplyTransIR ->
            template <class MIP>
            static void ApplyTransAdd (const FiniteElement & fel,
                    const MIP & mip,
                    FlatVector<Complex> x, BareSliceVector<Complex> y,
                    LocalHeap & lh)
            {
                if(fel.ComplexShapes())
                {
                cout << __FILE__<<" "<<__LINE__<<endl;
      FlatMatrixFixHeight<DIM_DMAT, Complex> mat(DIM*fel.GetNDof(), lh);
      y.Range(DIM*fel.GetNDof()) += (Cast(fel).GetDShapeComplex(mip,lh)) * x;
      //GenerateMatrix (fel, mip, mat, lh);
      //y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                }
                else
                {
      FlatMatrixFixHeight<DIM_DMAT, double> mat(DIM*fel.GetNDof(), lh);
      GenerateMatrix (fel, mip, mat, lh);
      y.Range(DIM*fel.GetNDof()) += Trans (mat) * x;
                }
            }

    template <typename MIR, class TVX, class TVY>
    static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
                  const TVX & x, TVY & y,
                  LocalHeap & lh)
    {
                cout << __FILE__<<" "<<__LINE__<<endl;
      y.Range(0,DIM*fel.GetNDof()) = 0.0;
      for (size_t i = 0; i < mir.Size(); i++)
        {
          HeapReset hr(lh);
          ApplyTransAdd (fel, mir[i], x.Row(i), y, lh);
        }
    }

            //template <class MIR>
            //static void ApplyTransIR (const FiniteElement & fel,
                    //const MIR & mir,
                    //FlatMatrix<Complex> x, BareSliceVector<Complex> y,
                    //LocalHeap & lh)
            //{
                //cout << __FILE__<<" "<<__LINE__<<endl;
                ////DiffOp<DiffOpMappedGradient<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);
            //}

            using DiffOp<DiffOpMappedGradient<D, FEL> >::AddTransSIMDIR;
            static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
            {
                cout << __FILE__<<" "<<__LINE__<<endl;
                Cast(fel).AddGradTrans (mir, y, x);
            }

    };


  template <int D>
  class DiffOpMappedHesse : public DiffOp<DiffOpMappedHesse<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 2 };

    static string Name() { return "hesse"; }
    // static Array<int> GetDimensions() { return Array<int> ( { D,D } ); }
    static INT<2> GetDimensions() { return { D,D }; }

    static auto & Cast (const FiniteElement & fel)
    { return static_cast<const ScalarMappedElement<D>&> (fel); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      Cast(fel).CalcMappedDDShape(mip, Trans(mat));
    }
  };


}
#endif
