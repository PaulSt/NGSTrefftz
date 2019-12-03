#ifndef FILE_TREFFTZGPPWELEMENT_HPP
#define FILE_TREFFTZGPPWELEMENT_HPP

#include <fem.hpp>
#include "helpers.hpp"
#include "scalarmappedfe.hpp"

namespace ngfem
{
    template <int D>
        class TrefftzGppwFE : public ScalarMappedElement<D+1>
    {
        private:
            const int ord;
            const int npoly;
            Vec<D+1> elcenter;
            double elsize;
            float c;
            ELEMENT_TYPE eltype;
            int basistype;
            int gppword;
            const Array<double> &gamma;

        public:
            TrefftzGppwFE(const Array<double> &agamma, int agppword = 0, int aord = 1, float ac = 1.0, Vec<D+1> aelcenter = 0, double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG, int abasistype = 0);

            float GetWavespeed() const { return c; }
            void SetWavespeed(double wavespeed) {c = wavespeed;}

            //TrefftzWaveFE<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return this;}
            //TrefftzWaveFE<D> * SetElSize(double aelsize) {elsize = aelsize; return this;}
            // TrefftzWaveFE<D> * SetWavespeed(float ac) {c = ac; return this;}

            virtual ELEMENT_TYPE ElementType() const { return eltype; }

            using ScalarMappedElement<D+1>::CalcShape;
            virtual void CalcShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const;
            virtual void CalcShape (const SIMD_BaseMappedIntegrationRule & smir, BareSliceMatrix<SIMD<double>> shape) const;

            using ScalarMappedElement<D+1>::CalcDShape;
            virtual void CalcDShape (const BaseMappedIntegrationPoint & mip, BareSliceMatrix<> dshape) const;
            virtual void CalcDShape (const SIMD_BaseMappedIntegrationRule & smir, BareSliceMatrix<SIMD<double>> dshape) const;

            int GetNDof() const {return this->ndof;}

            static int NDirections(int ord) { return ord==0?1:2;}
            static double GetDirection(int ord, int k) { return k%2?1:-1;}

            //using ScalarMappedElement<D>::CalcMappedDShape;
            //HD NGS_DLL_HEADER virtual
            void CalcMappedDShape (const BaseMappedIntegrationPoint & bmip,
                    BareSliceMatrix<> dshape) const
            { CalcDShape(bmip,dshape); }
            void CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> dshapes) const
            { CalcDShape(mir,dshapes); }

            void Evaluate (const SIMD_BaseMappedIntegrationRule & mir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const
            {
                STACK_ARRAY(SIMD<double>, mem, this->ndof*mir.Size());
                FlatMatrix<SIMD<double>> shape(this->ndof,mir.Size(),&mem[0]);
                CalcShape (mir, shape);
                const int nsimd = SIMD<double>::Size();
                FlatMatrix<double> bdbmat(this->ndof,mir.Size()*nsimd,&shape(0,0)[0]);
                FlatVector<double> bdbvec(mir.Size()*nsimd,&values(0)[0]);
                bdbvec = Trans(bdbmat) * coefs;
            }
            void AddTrans (const SIMD_BaseMappedIntegrationRule & mir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const
            {
                STACK_ARRAY(SIMD<double>, mem, this->ndof*mir.Size());
                FlatMatrix<SIMD<double>> shape(this->ndof,mir.Size(),&mem[0]);
                CalcShape (mir, shape);
                const int nsimd = SIMD<double>::Size();
                FlatMatrix<double> bdbmat(this->ndof,mir.Size()*nsimd,&shape(0,0)[0]);
                FlatVector<double> bdbvec(mir.Size()*nsimd,&values(0)[0]);
                coefs.AddSize(this->ndof) += bdbmat * bdbvec;
            }

            void EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
            {
                STACK_ARRAY(SIMD<double>, mem, (D+1)*this->ndof*ir.Size());
                FlatMatrix<SIMD<double>> simddshapes((D+1)*this->ndof,ir.Size(),&mem[0]);
                CalcDShape(ir,simddshapes);
                const int nsimd = SIMD<double>::Size();
                FlatMatrix<double> dshapes(this->ndof,(D+1)*nsimd*ir.Size(),&simddshapes(0,0)[0]);
                FlatVector<double> bdbvec((D+1)*nsimd*ir.Size(),&values(0,0)[0]);
                bdbvec = Trans(dshapes)*coefs;
            }
            void AddGradTrans (const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> values,
                    BareSliceVector<> coefs) const
            {
                STACK_ARRAY(SIMD<double>, mem, (D+1)*this->ndof*mir.Size());
                FlatMatrix<SIMD<double>> simddshapes((D+1)*this->ndof,mir.Size(),&mem[0]);
                CalcDShape(mir,simddshapes);
                const int nsimd = SIMD<double>::Size();
                FlatMatrix<double> dshapes(this->ndof,(D+1)*nsimd*mir.Size(),&simddshapes(0,0)[0]);
                FlatVector<double> bdbvec((D+1)*nsimd*mir.Size(),&values(0,0)[0]);
                coefs.AddSize(this->ndof) += dshapes*bdbvec;
            }

    };


    template<int D>
    class TrefftzGppwBasis{
        public:
            static TrefftzGppwBasis& getInstance(){
                static TrefftzGppwBasis ginstance;
                // volatile int dummy{};
                return ginstance;
            }

            const CSR* TB(int ord, int gppword, const Array<double> &gamma, int basistype = 0);

        private:
            TrefftzGppwBasis()= default;
            ~TrefftzGppwBasis()= default;
            TrefftzGppwBasis(const TrefftzGppwBasis&)= delete;
            TrefftzGppwBasis& operator=(const TrefftzGppwBasis&)= delete;

            mutex gentrefftzbasis;
            std::map<std::string,CSR> gtbstore;
            //once_flag tbonceflag;
            //void TB_inner(const Array<double> &gamma, int ord, Matrix<> &trefftzbasis, Vec<D+1, int> coeffnum, int basis, int dim, int &tracker, int basistype);
            void TB_inner(int ord, Matrix<> &trefftzbasis, Vec<D+1, int> coeffnum, int basis, int dim, int &tracker, int basistype, double wavespeed);
            int IndexMap2(Vec<D+1, int> index, int ord);
    };

}


//#ifdef NGS_PYTHON
//#include <python_ngstd.hpp>
//void ExportTrefftzElement(py::module m);
//#endif // NGS_PYTHON


#endif // FILE_TrefftzGPPWElement_HPP