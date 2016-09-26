#include "base.h"
#include "gtest/gtest.h"
using namespace std;
using testing::Types;

typedef Types<agl::xorshift64star, agl::xorshift1024star> RandomTestTypes;

template<typename T>
class random_test : public testing::Test {};
TYPED_TEST_CASE(random_test, RandomTestTypes);

TYPED_TEST(random_test, pi) {
  constexpr size_t kNumPoints = 1000000;
  TypeParam rng;

  size_t c = 0;
  for (size_t i = 0; i < kNumPoints; ++i) {
    uniform_real_distribution<double> urd(0.0, 1.0);
    const double x = urd(rng);
    const double y = urd(rng);
    if (x * x + y * y < 1.0) ++c;
  }

  const double pi = c * 4.0 / kNumPoints;
  ASSERT_LE(3.14, pi);
  ASSERT_LE(pi, 3.15);
}
