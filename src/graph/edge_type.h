#pragma once
#include <iostream>
#include "graph.h"

namespace agl {
//
// Unweighted edge
//

// By directly making |unweighted_edge| an alias to |V|
// (instead of making a new struct),
// edge lists for unweighted graphs become |std::pair<V, V>|,
// which is useful and beautiful to handle them.
using unweighted_edge = V;

constexpr V to(const unweighted_edge &e) {
  return e;
}

// The return type of this function determines the type of |W|.
constexpr int32_t weight(const unweighted_edge&) {
  return 1;
}

constexpr unweighted_edge reverse_edge(V v, const unweighted_edge&) {
  return unweighted_edge{v};
}

//
// Weighted edge
//
template<typename WeightType>
struct weighted_edge {
  using weight_type = WeightType;

  V to;
  weight_type weight;
};

template<typename WeightType>
constexpr V to(const weighted_edge<WeightType> &e) {
  return e.to;
}

template<typename WeightType>
constexpr WeightType weight(const weighted_edge<WeightType> &e) {
  return e.weight;
}

template<typename WeightType>
constexpr weighted_edge<WeightType> reverse_edge(V v, const weighted_edge<WeightType> &e) {
  return weighted_edge<WeightType>{v, e.weight};
}

template<typename WeightType>
std::ostream &operator<<(std::ostream &os, const weighted_edge<WeightType> &e) {
  return os << e.to << "(" << e.weight << ")";
};

template<typename WeightType>
std::istream &operator>>(std::istream &is, weighted_edge<WeightType> &e) {
  return is >> e.to >> e.weight;
};
}  // namespace agl
