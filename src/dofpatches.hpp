#ifndef FILE_DOFPATCHES_HPP
#define FILE_DOFPATCHES_HPP

#include <comp.hpp>

namespace ngcomp
{
  /// Independent trial DOF groups; their element supports may overlap.
  /// Matching trial/test providers must number corresponding groups alike.
  struct DofPatches
  {
    Table<size_t> elements;
    Table<DofId> dofs;
  };

  class DofPatchProvider
  {
  public:
    virtual ~DofPatchProvider () = default;
    virtual DofPatches GetDofPatches () const = 0;
  };
}

#endif
