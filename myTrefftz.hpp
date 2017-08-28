#ifndef FILE_MYTREFFTZ_HPP
#define FILE_MYTREFFTZ_HPP

#include "MultiArray.hpp"

namespace ngfem
{
  template <int D> class MyTrefftz : public FiniteElement
  {
  private:
    const int order;
    vector<MultiArray<float, D + 1>> basisFunctions;
    int nbasis;

  public:
    // constructor
    MyTrefftz () : FiniteElement () { ; }

    MyTrefftz (int aorder) : order (aorder)
    {
      nbasis
          = BinCoeff (D + order, order) + BinCoeff (D + order - 1, order - 1);
      basisFunctions = vector<MultiArray<float, D + 1>> (
          nbasis, aorder + 1); // nbasis MultiArrays of depth aorder+1
      cout << "order: " + to_string (order) + ", dimension: " + to_string (D)
                  + "\n";
    }

    virtual ELEMENT_TYPE ElementType () const { return ET_TRIG; }

    virtual void
    CalcShape (const IntegrationPoint &ip, BareSliceVector<> shape) const;

    virtual void
    CalcDShape (const IntegrationPoint &ip, SliceMatrix<> dshape) const;

    void TrefftzBasis ();

    void MakeIndices (int maxes, vector<array<int, D + 1>> &indices);

    void MakeIndices_inner (int dim, array<int, D + 1> &numbers, int maxes,
                            vector<array<int, D + 1>> &result);

    int BinCoeff (int n, int k) const;

    constexpr int ipow_ar (array<int, D + 1> base, array<int, D + 1> exp)
    {
      int result = 1;
      for (int i = 0; i < D + 1; i++)
        {
          result *= ipow (base[i], exp[i]);
        }
      return result;
    }
  };
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMyTrefftz (py::module m);
#endif // NGS_PYTHON

#endif // FILE_MYTREFFTZ_HPP
