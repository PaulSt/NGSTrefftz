#ifndef FILE_MULTIARRAY_HPP
#define FILE_MULTIARRAY_HPP

//#include <iostream>
//#include <cmath>
//#include <array>
//#include <cstdint>
// using namespace std;
#include <vector>

#include "helpers.cpp"

namespace ngfem
{
  template <class T, int depth> class MultiArray
  {
  private:
    std::vector<T> multiArray;
    int each_size;

  public:
    MultiArray () { each_size = 0; }

    MultiArray (int aeach_size) //: each_size(aeach_size)
    {
      this->each_size = aeach_size;
      // cout << "MultiArray initilized with size: " << ipow(aeach_size, depth)
      // << " dim: " << depth << " order: " << aeach_size << endl;
      multiArray.resize (ipow (aeach_size, depth));
    }

    void put (array<int, depth> indices, T thing_to_put);

    T get (array<int, depth> indices) const;

    T operator[] (const std::initializer_list<int> indices);
  };

}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMultiArray (py::module m);
#endif // NGS_PYTHON

#endif // FILE_MULTIARRAY_HPP
