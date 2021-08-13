#ifndef FILE_TREFFTZWAVEELEMENT_HPP
#define FILE_TREFFTZWAVEELEMENT_HPP

#include <fem.hpp>
#include "scalarmappedfe.hpp"

namespace ngfem
{
  template <int D> class TrefftzWaveFE : public ScalarMappedElement<D + 1>
  {
  private:
    const int ord;
    const int npoly;
    ELEMENT_TYPE eltype;
    int basistype;

  public:
    // TrefftzWaveFE();
    TrefftzWaveFE (int aord = 1, float ac = 1.0, Vec<D + 1> aelcenter = 0,
                   double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG,
                   int abasistype = 0);

    double GetWavespeed () const { return this->c; }
    void SetWavespeed (double wavespeed) { this->c = wavespeed; }

    // TrefftzWaveFE<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return
    // this;} TrefftzWaveFE<D> * SetElSize(double aelsize) {elsize = aelsize;
    // return this;}
    //  TrefftzWaveFE<D> * SetWavespeed(float ac) {c = ac; return this;}
    //
    ELEMENT_TYPE ElementType () const { return eltype; }

    using ScalarMappedElement<D + 1>::CalcShape;
    using ScalarMappedElement<D + 1>::CalcDShape;

    using ScalarMappedElement<D + 1>::CalcMappedDShape;

    void CalcMappedDDShape (const BaseMappedIntegrationPoint &bmip,
                            BareSliceMatrix<> hddshape) const
    {
      cout << "hesse not implemented" << endl;
    }
  };

  template <int D> class TrefftzWaveBasis
  {
  public:
    static TrefftzWaveBasis &getInstance ()
    {
      static TrefftzWaveBasis instance;
      // volatile int dummy{};
      return instance;
    }

    const CSR TB (int ord);
    void CreateTB (int ord, int basistype = 0);
    static void TB_inner (int ord, Matrix<> &trefftzbasis,
                          Vec<D + 1, int> coeffnum, int basis, int dim,
                          int &tracker, int basistype, double wavespeed = 1.0);
    static int IndexMap2 (Vec<D + 1, int> index, int ord);

  private:
    TrefftzWaveBasis () = default;
    ~TrefftzWaveBasis () = default;
    TrefftzWaveBasis (const TrefftzWaveBasis &) = delete;
    TrefftzWaveBasis &operator= (const TrefftzWaveBasis &) = delete;

    Array<CSR> tbstore;
    // once_flag tbonceflag;
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement (py::module m);
#endif // NGS_PYTHON

#endif // FILE_TrefftzWaveElement_HPP
