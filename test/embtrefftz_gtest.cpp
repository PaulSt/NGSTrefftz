#include <gtest/gtest.h>

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

TEST (max_a_minus_b_or_0_Test, aBiggerThanB)
{
  ASSERT_EQ (max_a_minus_b_or_0 (3, 2), 1);
}

TEST (max_a_minus_b_or_0_Test, aEqualsB)
{
  ASSERT_EQ (max_a_minus_b_or_0 (4, 4), 0);
}

TEST (max_a_minus_b_or_0_Test, aLessThanBNoUnderflow)
{
  ASSERT_EQ (max_a_minus_b_or_0 (3, 4), 0);
}

TEST (calcNdofTrefftzTest, noTrefftzOp)
{
  ASSERT_EQ (
      calcNdofTrefftz (3, 3, 1, size_t (3), true, SliceVector<double>{}), 2);
}

TEST (calcNdofTrefftzTest, noTrefftzOpUnderflow)
{
  ASSERT_EQ (
      calcNdofTrefftz (3, 3, 4, size_t (3), true, SliceVector<double>{}), 0);
}

TEST (calcNdofTrefftzTest, NdofTrefftzGiven)
{
  ASSERT_EQ (
      calcNdofTrefftz (3, 4, 5, size_t (42), false, SliceVector<double>{}),
      42);
}

TEST (calcNdofTrefftzTest, EpsilonGiven)
{
  ASSERT_EQ (
      calcNdofTrefftz (42, 4, 5, double (1e-2), false, SliceVector<double>{}),
      42 - 4 - 5);
}

TEST (calcNdofTrefftzTest, EpsilonGivenWithSingularValues)
{
  // two extra trefftz dofs from the singuar vlaues
  // with eps == 1e-2.
  auto singular_values = Vector<double> (3);
  singular_values[0] = 1e2;
  singular_values[1] = 1e-4;
  singular_values[2] = 1e-3;

  ASSERT_EQ (calcNdofTrefftz<double> (42, 4, 5, double (1e-2), false,
                                      singular_values.View ()),
             42 - 4 - 5 + 2);
}

TEST (copyBitArrayTest, CopySourceZeros)
{
  auto source = make_shared<BitArray> (2);
  auto target = make_shared<BitArray> (10);

  source->Clear ();
  target->Set ();

  copyBitArray (target, source);

  for (size_t i = 0; i < source->Size (); i++)
    EXPECT_EQ (target->Test (i), false);
  for (size_t i = source->Size (); i < target->Size (); i++)
    EXPECT_EQ (target->Test (i), true);
}

TEST (copyBitArrayTest, CopySourceOnes)
{
  auto source = make_shared<BitArray> (3);
  auto target = make_shared<BitArray> (10);

  source->Set ();
  target->Clear ();

  copyBitArray (target, source);

  for (size_t i = 0; i < source->Size (); i++)
    EXPECT_EQ (target->Test (i), true);
  for (size_t i = source->Size (); i < target->Size (); i++)
    EXPECT_EQ (target->Test (i), false);
}
