#pragma once
// https://ideone.com/rj3tDd
// https://twitter.com/ir5/status/435686310792007680

namespace agl {
template<typename T>
class irange {
public:
  struct iterator_type {
    T x;
    T operator*() const { return x; }
    bool operator!=(const iterator_type& lhs) const { return x < lhs.x; }
    void operator++() { ++x; }
  };

  irange(T n) : i_(), n_{n} {}
  iterator_type begin() { return i_; }
  iterator_type end() { return n_; }

private:
  iterator_type i_, n_;
};

template<typename T>
irange<T> make_irange(T n) {
  return irange<T>(n);
}
}  // namespace agl
