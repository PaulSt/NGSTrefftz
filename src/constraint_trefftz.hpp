#include <expr.hpp>
#include <fespace.hpp>
#include <memory>
#include <fem.hpp>

namespace ngcomp
{
  template <class SCLA>
  std::shared_ptr<ngbla::Matrix<SCLA>>
  ConstraintTrefftz (std::shared_ptr<ngfem::SumOfIntegrals> op,
                     std::shared_ptr<FESpace> fes,
                     std::shared_ptr<ngfem::SumOfIntegrals> cop_lhs,
                     std::shared_ptr<ngfem::SumOfIntegrals> cop_rhs,
                     std::shared_ptr<FESpace> fes_constraint, int trefftzndof);
}
