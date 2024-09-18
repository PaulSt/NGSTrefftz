#include <fespace.hpp>
#include <gtest/gtest.h>

#include <comp.hpp>
#include <stdexcept>
#include "embtrefftz.cpp"

TEST (reorderMatrixColumnsTest, throwOnDimensionMismatch)
{
  LocalHeap lh (1000);

  FlatMatrix<double> mat (1, 2, lh);
  Array<DofId> empty_arr (1);

  ASSERT_THROW (reorderMatrixColumns (mat, empty_arr, lh),
                std::invalid_argument);
}

TEST (reorderMatrixColumnsTest, doesNotFailOnEmptyInput)
{
  LocalHeap lh (1000);

  FlatMatrix<double> mat (0, 0, lh);
  Array<DofId> empty_arr (0);

  reorderMatrixColumns (mat, empty_arr, lh);
}

TEST (reorderMatrixColumnsTest, workOnSingularInput)
{
  LocalHeap lh (1000);

  FlatMatrix<int> mat (1, 1, lh);
  FlatMatrix<int> mat_expected (1, 1, lh);
  Array<DofId> arr (1);

  arr[0] = 1;
  mat (0, 0) = 3;
  mat_expected (0, 0) = 3;

  reorderMatrixColumns (mat, arr, lh);
  ASSERT_EQ (mat (0, 0), mat_expected (0, 0));
}

TEST (reorderMatrixColumnsTest, sortProperly)
{
  LocalHeap lh (1000);

  FlatMatrix<int> mat (1, 3, lh);
  FlatMatrix<int> mat_expected (1, 3, lh);
  Array<DofId> arr (3);

  arr[0] = 3;
  arr[1] = 2;
  arr[2] = 1;
  mat (0, 0) = 4;
  mat (0, 1) = 5;
  mat (0, 2) = 6;
  mat_expected (0, 0) = 6;
  mat_expected (0, 1) = 5;
  mat_expected (0, 2) = 4;

  reorderMatrixColumns (mat, arr, lh);
  for (int j : { 0, 1, 2 })
    {
      EXPECT_EQ (mat (0, j), mat_expected (0, j));
    }
}
