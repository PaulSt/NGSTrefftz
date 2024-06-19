#ifndef FILE_TREFFTZ_HELPER_HPP
#define FILE_TREFFTZ_HELPER_HPP

#include <core/localheap.hpp>
#include <elementtopology.hpp>
#include <expr.hpp>
#include <fem.hpp>
#include <fespace.hpp>
#include <iostream>
#include <matrix.hpp>
#include <memory>
#include <meshaccess.hpp>
#include <meshing/localh.hpp>
#include <meshing/meshing3.hpp>
#include <ostream>
#include <sparsematrix.hpp>
#include <vector>

using namespace ngfem;
using namespace ngcomp;

/// @param bf symblic representation of a bilinear form
/// @param bfis stores the calculated BilinearFormIntegrators
void calculateBilinearFormIntegrators (
    const SumOfIntegrals &bf,
    Array<shared_ptr<BilinearFormIntegrator>> bfis[4]);

/// @param lf symblic representation of a linear form
/// @param lfis stores the calculated LinearFormIntegrators
void calculateLinearFormIntegrators (
    SumOfIntegrals &lf, Array<shared_ptr<LinearFormIntegrator>> lfis[4]);

/// decides, if the given finite element space
/// has hidden degrees of freedom.
///
/// @param fes finite element space
bool fesHasHiddenDofs (const FESpace &fes);

/// Tests, if the given bilinear form is defined on the given element.
///
/// @param bf bilinear form
/// @param mesh_element local mesh element
bool bfIsDefinedOnElement (const SumOfIntegrals &bf,
                           const Ngs_Element &mesh_element);

/// adds the integration of `bf_integrators` to the element local matrix of a
/// bilinear form.
///
/// @param elmat element matrix is written into this matrix. Should be of
///     dimension (ndof, ndof), where ndof is the number of local dofs of the
///     finite element space on this element.
/// @param bf_integrators array of integrators, that correspond to a bilinear
///     form. Here, only the volume integrators are relevant.
/// @param mesh_access access to the mesh Elements
/// @param element_id id of the local element
/// @param fes trial finite element space of the bilinear form
/// @param test_fes test finite element space of the bilinear form. May be the
///     same as `fes`.
/// @param local_heap allocator. elmat will be allocated here, as well as
template <class SCAL>
inline void addIntegrationToElementMatrix (
    FlatMatrix<SCAL> elmat,
    const Array<shared_ptr<BilinearFormIntegrator>> &bf_integrators,
    const MeshAccess &mesh_access, const ElementId &element_id,
    const FESpace &fes, const FESpace &test_fes, LocalHeap &local_heap)
{
  const HeapReset hr (local_heap);

  auto &trafo = mesh_access.GetTrafo (element_id, local_heap);

  auto &test_fel = test_fes.GetFE (element_id, local_heap);
  auto &trial_fel = fes.GetFE (element_id, local_heap);

  const bool mixed_mode = std::addressof (test_fes) != std::addressof (fes);

  bool symmetric_so_far = true;
  int bfi_ind = 0;
  while (bfi_ind < bf_integrators.Size ())
    {
      auto &bfi = bf_integrators[bfi_ind];
      bfi_ind++;
      if (bfi->DefinedOnElement (element_id.Nr ()))
        {
          auto &mapped_trafo = trafo.AddDeformation (
              bfi->GetDeformation ().get (), local_heap);
          try
            {
              if (mixed_mode)
                {
                  const auto &mixed_fel
                      = MixedFiniteElement (trial_fel, test_fel);
                  bfi->CalcElementMatrixAdd (mixed_fel, mapped_trafo, elmat,
                                             symmetric_so_far, local_heap);
                }
              else
                {
                  bfi->CalcElementMatrixAdd (test_fel, mapped_trafo, elmat,
                                             symmetric_so_far, local_heap);
                }
            }
          catch (ExceptionNOSIMD e)
            {
              elmat = 0.0;
              cout << IM (6) << "ExceptionNOSIMD " << e.What () << endl
                   << "switching to scalar evaluation" << endl;
              bfi->SetSimdEvaluate (false);
              bfi_ind = 0;
            }
        }
    }
}

/// extract only the visible dofs
/// @param elmat assembled element matrix of a bilinear form
/// @param element_id id of the local element
/// @param fes trial fe space of the bilinear form
/// @param fes test fe space of the bilinear form
/// @param dofs dofs of `fes`
/// @param test_dofs dofs of `test_fes`
/// @param local_heap allocator
/// @tparam SCAL scalar value type and matrix dimensions
template <class SCAL>
void extractVisibleDofs (FlatMatrix<SCAL> elmat, const ElementId &element_id,
                         const FESpace &fes, const FESpace &test_fes,
                         Array<DofId> &dofs, Array<DofId> &test_dofs,
                         LocalHeap &local_heap)
{
  const HeapReset hr (local_heap);

  Array<DofId> vdofs, vtest_dofs;
  fes.GetDofNrs (element_id, vdofs, VISIBLE_DOF);
  test_fes.GetDofNrs (element_id, vtest_dofs, VISIBLE_DOF);
  FlatMatrix<SCAL> velmat (vtest_dofs.Size (), vdofs.Size (), local_heap);
  for (int jj = 0; jj < dofs.Size (); jj++)
    for (int ii = 0; ii < test_dofs.Size (); ii++)
      {
        auto j = vdofs.Pos (dofs[jj]);
        auto i = vtest_dofs.Pos (test_dofs[ii]);
        if (i != size_t (-1) && j != size_t (-1))
          velmat (i, j) = elmat (ii, jj);
      }
  dofs = std::move (vdofs);
  test_dofs = std::move (vtest_dofs);
  elmat.Assign (velmat);
}
#endif
