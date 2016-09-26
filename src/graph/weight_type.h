#pragma once
#include <algorithm>
#include <random>
#include "base/base.h"

namespace agl {
//
// W traits
//
template<class T>
struct weight_traits {
  static constexpr inline T eps() {
    return T();
  }

  static constexpr inline T abs(T x) {
    return std::abs(x);
  }

  static constexpr inline bool is_zero(T x) {
    return x == T();
  }

  static constexpr inline bool is_equal(T x, T y) {
    return x == y;
  }

  static constexpr inline T random() {
    return std::uniform_int_distribution<T>(0, 10)(agl::random);
  }

  static constexpr inline T unit() {
    return 1;
  }
};

template<class T, class Derived>
struct weight_traits_decimal {
  static constexpr inline T abs(T x) {
    return std::abs(x);
  }

  static constexpr inline bool is_zero(T x) {
    return abs(x) <= Derived::eps();
  }

  static constexpr inline bool is_equal(T x, T y) {
    return is_zero(y - x) || is_zero((y - x) / std::max(abs(x), abs(y)));
  }

  static constexpr inline T random() {
    return std::uniform_real_distribution<T>(0.0, 1.0)(agl::random);
  }

  static constexpr inline T unit() {
    return 1;
  }
};

template<> struct weight_traits<float>
        : public weight_traits_decimal<float, weight_traits<float>> {
  static constexpr inline float eps() {
    return 1E-4F;
  }
};

template<> struct weight_traits<double>
        : public weight_traits_decimal<double, weight_traits<double>> {
  static constexpr inline double eps() {
    return 1E-9;
  }
};

template<> struct weight_traits<long double>
        : public weight_traits_decimal<long double, weight_traits<long double>> {
  static constexpr inline long double eps() {
    return 1E-12L;
  }
};

//
// Interfaces
//
template<typename T>
constexpr inline T abs(T x) {
  return weight_traits<T>::abs(x);
}

template<typename T>
constexpr inline bool is_zero(T x) {
  return weight_traits<T>::is_zero(x);
}

template<typename T>
constexpr inline int signum(T x) {
  return is_zero(x) ? T() : (x > T() ? +1 : -1);
}

template<typename T>
inline bool is_eq(T x, T y) {
  return weight_traits<T>::is_equal(x, y);
}

template<typename T>
inline bool is_lt(T x, T y) {
  return is_eq(x, y) ? false : (x < y);
}

template<typename T>
inline bool is_le(T x, T y) {
  return !is_lt(y, x);
}

template<typename T>
inline T random_weight() {
  return weight_traits<T>::random();
}

template<typename T>
inline T unit_weight() {
  return weight_traits<T>::unit();
}
}  // namespace agl
