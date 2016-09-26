#include "weight_type.h"
#include "gtest/gtest.h"
using namespace agl;

TEST(WeightType, Trivial) {
  EXPECT_TRUE(is_zero(0));
  EXPECT_FALSE(is_zero(1));

  // float
  EXPECT_TRUE(is_zero(+0.0f));
  EXPECT_TRUE(is_zero(+1E-11f));
  EXPECT_TRUE(is_zero(-1E-11f));
  EXPECT_TRUE(is_zero(+1E-6f));
  EXPECT_TRUE(is_zero(-1E-6f));
  EXPECT_FALSE(is_zero(+1E-3f));
  EXPECT_FALSE(is_zero(-1E-3f));

  // double
  EXPECT_TRUE(is_zero(+0.0));
  EXPECT_TRUE(is_zero(+1E-11));
  EXPECT_TRUE(is_zero(-1E-11));
  EXPECT_FALSE(is_zero(+1E-6));
  EXPECT_FALSE(is_zero(-1E-6));
  EXPECT_FALSE(is_zero(+1E-3));
  EXPECT_FALSE(is_zero(-1E-3));

  // long double
  EXPECT_TRUE(is_zero(+0.0L));
  EXPECT_FALSE(is_zero(+1E-11L));
  EXPECT_FALSE(is_zero(-1E-11L));
  EXPECT_FALSE(is_zero(+1E-6L));
  EXPECT_FALSE(is_zero(-1E-6L));
  EXPECT_FALSE(is_zero(+1E-3L));
  EXPECT_FALSE(is_zero(-1E-3L));
}
