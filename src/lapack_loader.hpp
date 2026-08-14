#ifndef NGSTREFFTZ_LAPACK_LOADER_HPP
#define NGSTREFFTZ_LAPACK_LOADER_HPP

#include <bla.hpp>

namespace ngstrefftz
{
  using RealGesdd
      = int (*) (char *, ngbla::integer *, ngbla::integer *, double *,
                 ngbla::integer *, double *, double *, ngbla::integer *,
                 double *, ngbla::integer *, double *, ngbla::integer *,
                 ngbla::integer *, ngbla::integer *);
  using ComplexGesdd
      = int (*) (char *, ngbla::integer *, ngbla::integer *, ngbla::Complex *,
                 ngbla::integer *, double *, ngbla::Complex *,
                 ngbla::integer *, ngbla::Complex *, ngbla::integer *,
                 ngbla::Complex *, ngbla::integer *, double *,
                 ngbla::integer *, ngbla::integer *);

  struct LapackGesdd
  {
    RealGesdd realgesdd;
    ComplexGesdd complexgesdd;
  };

  /// Resolves the LAPACK SVD entry points from NGSolve's loaded backend.
  const LapackGesdd &GetLapackGesdd ();
}

#endif
