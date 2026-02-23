#include <gtest/gtest.h>
#include "embtrefftz.hpp"

using namespace ngcomp;

TEST (calcNdofTrefftzTest, noTrefftzOp)
{
  ASSERT_EQ (calcNdofTrefftz<double> (3, 3, 1, size_t (3), true,
                                      SliceVector<double>{}),
             2);
}

TEST (calcNdofTrefftzTest, noTrefftzOpUnderflow)
{
  ASSERT_EQ (calcNdofTrefftz<double> (3, 3, 4, size_t (3), true,
                                      SliceVector<double>{}),
             0);
}

TEST (calcNdofTrefftzTest, NdofTrefftzGiven)
{
  ASSERT_EQ (calcNdofTrefftz<double> (3, 4, 5, size_t (42), false,
                                      SliceVector<double>{}),
             42);
}

TEST (calcNdofTrefftzTest, EpsilonGiven)
{
  ASSERT_EQ (calcNdofTrefftz<double> (42, 4, 5, double (1e-2), false,
                                      SliceVector<double>{}),
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
    EXPECT_FALSE (target->Test (i));
  for (size_t i = source->Size (); i < target->Size (); i++)
    EXPECT_TRUE (target->Test (i));
}

TEST (copyBitArrayTest, CopySourceOnes)
{
  auto source = make_shared<BitArray> (3);
  auto target = make_shared<BitArray> (10);

  source->Set ();
  target->Clear ();

  copyBitArray (target, source);

  for (size_t i = 0; i < source->Size (); i++)
    EXPECT_TRUE (target->Test (i));
  for (size_t i = source->Size (); i < target->Size (); i++)
    EXPECT_FALSE (target->Test (i));
}

TEST (TrefftzEmbeddingCtorTest, NdofAndEpsModeSelection)
{
  using std::make_shared;
  using std::shared_ptr;

  // Case 1: conflicting inputs => must throw invalid_argument immediately
  EXPECT_THROW (
      (ngcomp::TrefftzEmbedding (nullptr, nullptr, nullptr, nullptr,
                                 /*_ndof_trefftz=*/5,
                                 /*_eps=*/1e-8, nullptr, nullptr, nullptr,
                                 nullptr, nullptr, nullptr)),
      std::invalid_argument);

  // Case 2: valid arguments, but no FESpace, fail later
  try
    {
      ngcomp::TrefftzEmbedding emb (
          nullptr, nullptr, nullptr, nullptr,
          /*_ndof_trefftz=*/std::numeric_limits<size_t>::max (),
          /*_eps=*/1e-8, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

      FAIL () << "Expected constructor to throw due to missing FESpace";
    }
  catch (const std::invalid_argument &)
    {
      FAIL () << "Wrong branch taken: constructor treated sentinel "
                 "ndof_trefftz as conflict";
    }
  catch (...)
    {
      SUCCEED (); // expected: later failure (e.g. no fespace found)
    }
}
