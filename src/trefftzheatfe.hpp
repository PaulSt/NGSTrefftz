#ifndef FILE_TREFFTZHEATELEMENT_HPP
#define FILE_TREFFTZHEATELEMENT_HPP

#include <fem.hpp>
#include "scalarmappedfe.hpp"

namespace ngfem
{
  template <int D> class TrefftzHeatFE : public ScalarMappedElement<D + 1>
  {
  private:
    ELEMENT_TYPE eltype;
    int basistype;

  public:
    // TrefftzHeatFE();
    TrefftzHeatFE (int aord = 1, float ac = 1.0, Vec<D + 1> aelcenter = 0,
                   double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG,
                   int abasistype = 0);

    // double GetWavespeed() const { return this->c; }
    // void SetWavespeed(double wavespeed) {this->c = wavespeed;}

    // TrefftzWaveFE<D> * SetCenter(Vec<D> acenter) {elcenter = acenter; return
    // this;} TrefftzWaveFE<D> * SetElSize(double aelsize) {elsize = aelsize;
    // return this;}
    //  TrefftzWaveFE<D> * SetWavespeed(float ac) {c = ac; return this;}
    //
    virtual ELEMENT_TYPE ElementType () const { return eltype; }

    using ScalarMappedElement<D + 1>::CalcShape;
    using ScalarMappedElement<D + 1>::CalcDShape;

    using ScalarMappedElement<D + 1>::CalcMappedDShape;
    using ScalarMappedElement<D + 1>::CalcMappedDDShape;

    // void CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip,
    // BareSliceMatrix<> hddshape) const
    //{cout << "hesse not implemented" << endl;}
  };

  template <int D> class TrefftzHeatBasis
  {
  public:
    static TrefftzHeatBasis &getInstance ()
    {
      static TrefftzHeatBasis instance;
      // volatile int dummy{};
      return instance;
    }

    const CSR TB (int ord);
    void CreateTB (int ord, int basistype = 0);
    static void TB_inner (int ord, Matrix<> &trefftzbasis,
                          Vec<D + 1, int> coeffnum, int basis, int dim,
                          int &tracker, int basistype, double Heatspeed = 1.0);
    static int IndexMap2 (Vec<D + 1, int> index, int ord);

  private:
    TrefftzHeatBasis () = default;
    ~TrefftzHeatBasis () = default;
    TrefftzHeatBasis (const TrefftzHeatBasis &) = delete;
    TrefftzHeatBasis &operator= (const TrefftzHeatBasis &) = delete;

    Array<CSR> tbstore;
    // once_flag tbonceflag;
  };

  template <int D> class TrefftzHeatTestFE : public ScalarMappedElement<D + 1>
  {
  private:
    ELEMENT_TYPE eltype;
    int basistype;

  public:
    TrefftzHeatTestFE (int aord = 1, float ac = 1.0, Vec<D + 1> aelcenter = 0,
                       double aelsize = 1, ELEMENT_TYPE aeltype = ET_TRIG,
                       int abasistype = 0);

    virtual ELEMENT_TYPE ElementType () const { return eltype; }

    using ScalarMappedElement<D + 1>::CalcShape;
    using ScalarMappedElement<D + 1>::CalcDShape;

    using ScalarMappedElement<D + 1>::CalcMappedDShape;
    using ScalarMappedElement<D + 1>::CalcMappedDDShape;

    // void CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip,
    // BareSliceMatrix<> hddshape) const
    //{cout << "hesse not implemented" << endl;}
  };

  template <int D> class TrefftzHeatTestBasis
  {
  public:
    static TrefftzHeatTestBasis &getInstance ()
    {
      static TrefftzHeatTestBasis instance;
      // volatile int dummy{};
      return instance;
    }

    const CSR TB (int ord);
    void CreateTB (int ord, int basistype = 0);
    static void TB_inner (int ord, Matrix<> &trefftzbasis,
                          Vec<D + 1, int> coeffnum, int basis, int dim,
                          int &tracker, int basistype, double Heatspeed = 1.0);
    static int IndexMap2 (Vec<D + 1, int> index, int ord);

  private:
    TrefftzHeatTestBasis () = default;
    ~TrefftzHeatTestBasis () = default;
    TrefftzHeatTestBasis (const TrefftzHeatTestBasis &) = delete;
    TrefftzHeatTestBasis &operator= (const TrefftzHeatTestBasis &) = delete;

    Array<CSR> tbstore;
    // once_flag tbonceflag;
  };
}

#endif // FILE_TrefftzHeatElement_HPP
