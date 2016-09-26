/*
 * xorshift random generators, which are compatible to std::random
 * http://cpplover.blogspot.jp/2013/03/blog-post_22.html
 * http://xorshift.di.unimi.it/xorshift64star.c
 * http://xorshift.di.unimi.it/xorshift1024star.c
 */

#include <random>
#include <array>
#include <limits>
#include <gflags/gflags.h>
#include "macros.h"

DECLARE_int64(random_seed);

namespace agl {
class xorshift64star {
 public:
  typedef uint64_t result_type;

  static constexpr uint64_t min() {
    return std::numeric_limits<uint64_t>::min();
  }

  static constexpr uint64_t max() {
    return std::numeric_limits<uint64_t>::max();
  }

  xorshift64star(uint64_t seed = FLAGS_random_seed) : x_(seed) {
    CHECK(seed != 0);
  }

  uint64_t operator()() {
    x_ ^= x_ >> 12;
    x_ ^= x_ << 25;
    x_ ^= x_ >> 27;
    return x_ * 2685821657736338717LL;
  }

  uint64_t operator()(uint64_t limit) {  // [0, limit)
    return operator()() % limit;
  }

 private:
  uint64_t x_;
};

class xorshift1024star {
 public:
  typedef uint64_t result_type;

  static constexpr uint64_t min() {
    return std::numeric_limits<uint64_t>::min();
  }

  static constexpr uint64_t max() {
    return std::numeric_limits<uint64_t>::max();
  }

  // Delay the initialization to wait |FLAGS_random_seed| to be set.
  xorshift1024star() : p_(-1) {}

  xorshift1024star(uint64_t seed) {
    srand(seed);
  }

  void srand(uint64_t seed) {
    CHECK(seed != 0);
    p_ = 0;
    xorshift64star x(seed);
    for (int i = 0; i < 16; ++i) s_[i] = x();
  }


  uint64_t operator()() {
    if (p_ == -1) {
      srand(FLAGS_random_seed);
    }

    uint64_t s0 = s_[p_];
    uint64_t s1 = s_[p_ = (p_ + 1) & 15];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    return (s_[ p_ ] = s0 ^ s1) * 1181783497276652981LL;
  }

  uint64_t operator()(uint64_t limit) {  // [0, limit)
    return operator()() % limit;
  }

 private:
  int p_;
  std::array<uint64_t, 16> s_;
};

using random_type = xorshift1024star;

namespace {
// here's a global random number generator for lazy guys
random_type random;
}  // namespace
}  // namespace base
