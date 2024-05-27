#include <fem.hpp>
#include <fespace.hpp>
#include <integratorcf.hpp>
#include <memory>
#include <meshaccess.hpp>

using namespace ngfem;
using namespace ngcomp;

/// @param bf symblic representation of a bilinear form
/// @param bfis stores the calculated BilinearFormIntegrators
void calculateBilinearFormIntegrators (
    const SumOfIntegrals &bf,
    Array<shared_ptr<BilinearFormIntegrator>> bfis[4])
{
  for (auto icf : bf.icfs)
    {
      DifferentialSymbol &dx = icf->dx;
      bfis[dx.vb] += icf->MakeBilinearFormIntegrator ();
    }
}

/// @param lf symblic representation of a linear form
/// @param lfis stores the calculated LinearFormIntegrators
void calculateLinearFormIntegrators (
    SumOfIntegrals &lf, Array<shared_ptr<LinearFormIntegrator>> lfis[4])
{
  for (auto icf : lf.icfs)
    {
      DifferentialSymbol &dx = icf->dx;
      lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
    }
}

/// decides, if the given finite element space
/// has hidden degrees of freedom.
///
/// @param fes finite element space
bool fesHasHiddenDofs (const FESpace &fes)
{
  const size_t ndof = fes.GetNDof ();
  bool has_hidden_dofs = false;
  for (DofId d = 0; d < ndof; d++)
    if (HIDDEN_DOF == fes.GetDofCouplingType (d))
      return true;
  return false;
}

/// Tests, if the given bilinear form is defined on the given element.
///
/// @param bf bilinear form
/// @param mesh_element local mesh element
bool bfIsDefinedOnElement (const SumOfIntegrals &bf,
                           const Ngs_Element &mesh_element)
{
  for (auto icf : bf.icfs)
    {
      if (icf->dx.vb == VOL)
        if ((!icf->dx.definedonelements)
            || (icf->dx.definedonelements->Test (mesh_element.Nr ())))
          return true;
    }
  return false;
}
