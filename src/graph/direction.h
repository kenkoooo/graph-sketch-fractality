#pragma once
#include <type_traits>

namespace agl {
// We don't use scoped enum to avoid writing |static_cast| every time (we're lazy)
enum D : uint8_t {
  kFwd = 0,
  kBwd = 1,
  kNumDirections = 2,
};

using direction_underlying_type = std::underlying_type<D>::type;

class direction_range {
 public:
  struct iterator_type {
    direction_underlying_type x;
    D operator*() const { return static_cast<D>(x); }
    bool operator!=(const iterator_type &i) const { return x < i.x; }
    void operator++() { ++x; }
  };

  direction_range() : i_{0}, n_{kNumDirections} {}

  iterator_type begin() { return i_; }
  iterator_type end() { return n_; }

 private:
  iterator_type i_, n_;
};

inline direction_range directions() {
  return direction_range();
}

inline D reverse_direction(D d) {
  return D(1 - d);
}
}  // namespace agl
