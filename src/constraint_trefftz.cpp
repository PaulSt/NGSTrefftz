#include "constraint_trefftz.hpp"
#include <bilinearform.hpp>
#include <expr.hpp>
#include <memory>

namespace ngcomp
{
  template <class SCLA>
  std::shared_ptr<ngbla::Matrix<SCLA>>
  ConstraintTrefftz (std::shared_ptr<ngfem::SumOfIntegrals> op,
                     std::shared_ptr<FESpace> fes,
                     std::shared_ptr<ngfem::SumOfIntegrals> cop_lhs,
                     std::shared_ptr<ngfem::SumOfIntegrals> cop_rhs,
                     std::shared_ptr<FESpace> fes_constraint, int trefftzndof)
  {
    return make_shared<Matrix<SCLA>> ();
  }
}
