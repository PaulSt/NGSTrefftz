#ifndef FILE_SVDTREFFTZ_HPP
#define FILE_SVDTREFFTZ_HPP
#include <comp.hpp>
#include <fem.hpp>
#include <integratorcf.hpp>
#include <variant>

namespace ngcomp
{
  template <class SCAL>
  std::tuple<shared_ptr<BaseMatrix>,shared_ptr<BaseVector>> EmbTrefftz (shared_ptr<SumOfIntegrals> bf,
                       shared_ptr<FESpace> fes, 
                       shared_ptr<SumOfIntegrals> lf,
                       double eps, shared_ptr<FESpace> fes_test, int tndof
                   );
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportEmbTrefftz(py::module m);
#endif // NGS_PYTHON

#endif
