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

TEST (conformingTrefftzReshapeMatricesTest, viewsAreAssembled)
{
  LocalHeap local_heap (1000 * sizeof (int));
  size_t ndof = 9, ndof_conforming = 3, ndof_test = 5;
  FlatMatrix<int> elmat_A (4, 6, local_heap);
  FlatMatrix<int> elmat_B (5, 6, local_heap);
  FlatMatrix<int> elmat_Cl, elmat_L, elmat_Cr;

  conformingTrefftzReshapeMatrices (elmat_A, elmat_B, elmat_Cl, elmat_L,
                                    elmat_Cr, ndof, ndof_test,
                                    ndof_conforming);
  EXPECT_EQ (elmat_Cl.Shape (), make_tuple (ndof_conforming, ndof));
  EXPECT_EQ (elmat_L.Shape (), make_tuple (ndof_test, ndof));
  EXPECT_EQ (elmat_Cr.Shape (), make_tuple (ndof_conforming, ndof_conforming));
}

TEST (conformingTrefftzReshapeMatricesTest, viewsWork)
{
  LocalHeap local_heap (1000 * sizeof (int));
  size_t ndof = 4, ndof_conforming = 2, ndof_test = 3;
  FlatMatrix<int> elmat_A (ndof_conforming + ndof_test, ndof, local_heap);
  FlatMatrix<int> elmat_B (ndof_conforming + ndof_test, ndof_conforming,
                           local_heap);
  auto [elmat_Cl, elmat_L] = elmat_A.SplitRows (ndof_conforming);
  auto elmat_Cr = elmat_A.Rows (ndof_conforming);

  conformingTrefftzReshapeMatrices (elmat_A, elmat_B, elmat_Cl, elmat_L,
                                    elmat_Cr, ndof, ndof_test,
                                    ndof_conforming);

  size_t count = 0;
  for (int i = 0; i < elmat_A.Height (); i++)
    for (int j = 0; j < elmat_A.Width (); j++)
      elmat_A (i, j) = count++;
  count = 0;
  for (int i = 0; i < elmat_B.Height (); i++)
    for (int j = 0; j < elmat_B.Width (); j++)
      elmat_B (i, j) = count++;

  count = 0;
  for (int i = 0; i < elmat_Cl.Height (); i++)
    for (int j = 0; j < elmat_Cl.Width (); j++)
      {
        EXPECT_EQ (elmat_Cl (i, j), count++);
      }
  for (int i = 0; i < elmat_L.Height (); i++)
    for (int j = 0; j < elmat_L.Width (); j++)
      EXPECT_EQ (elmat_L (i, j), count++);

  count = 0;
  for (int i = 0; i < elmat_Cr.Height (); i++)
    for (int j = 0; j < elmat_Cr.Width (); j++)
      EXPECT_EQ (elmat_Cr (i, j), count++);
}

TEST (conformingTrefftzReshapeMatricesTest, shrinkConformingNdofCorrectly)
{
  LocalHeap local_heap (1000 * sizeof (int));
  size_t ndof = 4, ndof_conforming = 2, ndof_test = 3;
  FlatMatrix<int> elmat_A (ndof_conforming + ndof_test, ndof, local_heap);
  FlatMatrix<int> elmat_B (ndof_conforming + ndof_test, ndof_conforming,
                           local_heap);
  auto [elmat_Cl, elmat_L] = elmat_A.SplitRows (ndof_conforming);
  auto elmat_Cr = elmat_A.Rows (ndof_conforming);

  size_t count = 0;
  for (int i = 0; i < elmat_A.Height (); i++)
    for (int j = 0; j < elmat_A.Width (); j++)
      elmat_A (i, j) = count++;
  count = 0;
  for (int i = 0; i < elmat_B.Height (); i++)
    for (int j = 0; j < elmat_B.Width (); j++)
      elmat_B (i, j) = count++;

  ndof_conforming = 1;

  conformingTrefftzReshapeMatrices (elmat_A, elmat_B, elmat_Cl, elmat_L,
                                    elmat_Cr, ndof, ndof_test,
                                    ndof_conforming);

  ASSERT_EQ (elmat_A.Shape (), make_tuple (ndof_test + ndof_conforming, ndof));
  ASSERT_EQ (elmat_B.Shape (),
             make_tuple (ndof_test + ndof_conforming, ndof_conforming));
  ASSERT_EQ (elmat_Cl.Shape (), make_tuple (ndof_conforming, ndof));
  ASSERT_EQ (elmat_L.Shape (), make_tuple (ndof_test, ndof));
  ASSERT_EQ (elmat_Cr.Shape (), make_tuple (ndof_conforming, ndof_conforming));

  count = 0;
  for (int i = 0; i < elmat_Cl.Height (); i++)
    for (int j = 0; j < elmat_Cl.Width (); j++)
      {
        EXPECT_EQ (elmat_Cl (i, j), count++);
      }
  for (int i = 0; i < elmat_L.Height (); i++)
    for (int j = 0; j < elmat_L.Width (); j++)
      EXPECT_EQ (elmat_L (i, j), count++);

  count = 0;
  for (int i = 0; i < elmat_Cr.Height (); i++)
    for (int j = 0; j < elmat_Cr.Width (); j++)
      EXPECT_EQ (elmat_Cr (i, j), count++);
}
