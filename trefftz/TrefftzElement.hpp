#ifndef FILE_TREFFTZELEMENT_HPP
#define FILE_TREFFTZELEMENT_HPP

#include <fem.hpp>
#include <l2hofefo.hpp>
#include "helpers.cpp"
// using namespace ngfem;

namespace ngfem
{
  template <int D, int order>
  class TrefftzElement : public ScalarFiniteElement<D>
  {
  private:
    // vector<MultiArray<float,D+1> > basisFunctions;

    const int nbasis
        = BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1);
    array<array<int, D + 1>, BinCoeff (D + 1 + order, order)> indices;
    Matrix<double> basis;

  public:
    // constructor
    TrefftzElement ()
        : ScalarFiniteElement<D> (),
          basis (nbasis, BinCoeff (D + 1 + order,
                                   order)) //, npol(BinCoeff(D + order, order))
    {
      // basisFunctions = vector<MultiArray<float,D+1> >(nbasis, order+1);
      // //nbasis MultiArrays of depth order+1
      cout << "order: " + to_string (order) + ", dimension: " + to_string (D)
                  + ", number of basis functions: "
           << nbasis << endl;

      // cout << "\n ===== exponentials: \n";
      // indices.reserve( BinCoeff(D+1+order,order) );

      array<int, D + 1> numbers;

      int count = 0;
      MakeIndices (numbers, count);
      ~count;
      cout << "constructing basis" << endl;
      basis = 0;
      TrefftzBasis ();
    }

    virtual ELEMENT_TYPE ElementType () const { return ET_TRIG; }

    virtual void
    CalcShape (const IntegrationPoint &ip, BareSliceVector<> shape) const;

    virtual void
    CalcDShape (const IntegrationPoint &ip, SliceMatrix<> dshape) const;

    virtual void CalcShape (const BaseMappedIntegrationPoint &mip,
                            Vector<double> &shape) const;

    virtual void CalcDShape (const BaseMappedIntegrationPoint &mip,
                             SliceMatrix<> dshape) const;

    void TrefftzBasis ();

    void MakeIndices (array<int, D + 1> &numbers, int &count,
                      int dim = D + 1); //(int dim, array<int, D+1> &numbers,
                                        //vector< array<int, D+1> > &indices);

    float ipow_ar (IntegrationPoint base, array<int, D + 1> ex,
                   float result = 1, int count = D + 1) const;
    double ipow_ar (FlatVector<double> base, array<int, D + 1> ex,
                    float result = 1, int count = D + 1) const;

    int GetNBasis () const;

    int IndexMap (array<int, D + 1> index);
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportTrefftzElement (py::module m);
#endif // NGS_PYTHON

#endif // FILE_TrefftzElement_HPP
