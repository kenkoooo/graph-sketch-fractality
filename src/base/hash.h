#pragma once
#include <functional>

namespace std {  // argh!
template<typename A, typename B>
struct hash<std::pair<A, B>> {
  size_t operator()(const std::pair<A, B> &p) const {
    return std::hash<A>()(p.first) * 1000000007 + std::hash<B>()(p.second);
  }
};
}  // namespace std
