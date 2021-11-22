#ifndef FILE_SVDTREFFTZ_HPP
#define FILE_SVDTREFFTZ_HPP
#include <comp.hpp>
#include <fem.hpp>
#include <integratorcf.hpp>
#include <variant>

namespace ngcomp
{
  shared_ptr<BaseMatrix> SVDTrefftz (shared_ptr<SumOfIntegrals> bf,
                   shared_ptr<FESpace> fes, double eps=10e-10);
}

#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportSVDTrefftz(py::module m);
#endif // NGS_PYTHON

#endif
