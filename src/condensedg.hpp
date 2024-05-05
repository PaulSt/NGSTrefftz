#ifndef FILE_CONDENSEDG_HPP
#define FILE_CONDENSEDG_HPP
#include <comp.hpp>
#include <python_comp.hpp>
#include <fem.hpp>
// #include <integratorcf.hpp>
// #include <variant>
#include <bla.hpp>

namespace ngcomp
{
  void GetSubMatrix (shared_ptr<BaseMatrix> mat, FlatArray<int> drow,
                     FlatArray<int> dcol, FlatMatrix<> out);

  shared_ptr<BaseMatrix>
  CondenseDG (shared_ptr<BaseMatrix> mat, shared_ptr<BaseVector> vec,
              shared_ptr<FESpace> fes);
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportCondenseDG (py::module m);
#endif // NGS_PYTHON

#endif
