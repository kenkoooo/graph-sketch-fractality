#pragma once
#include "graph/graph.h"
#include "graph/weight_type.h"

namespace agl {
unweighted_edge_list generate_path(V num_vertices);
unweighted_edge_list generate_erdos_renyi(V num_vertices, double avg_deg);
unweighted_edge_list generate_ba(V final_num, V initial_num);
unweighted_edge_list generate_uv_flower(V required_num, V u, V v);
unweighted_edge_list generate_shm(V required_num, V initial_num, int t, double P = 0.0);

unweighted_edge_list make_undirected(const unweighted_edge_list &es);

template<typename GraphType>
typename GraphType::edge_list_type add_random_weight(const unweighted_edge_list &es);

template<>
unweighted_edge_list add_random_weight<unweighted_graph>(const unweighted_edge_list &es);

template<typename GraphType>
typename GraphType::edge_list_type add_random_weight(const unweighted_edge_list &es) {
  typename GraphType::edge_list_type wes(es.size());
  for (size_t i = 0; i < es.size(); ++i) {
    wes[i].first = es[i].first;
    wes[i].second = typename GraphType::E{es[i].second, random_weight<typename GraphType::W>()};
  }
  return wes;
}

template<typename GraphType>
typename GraphType::edge_list_type add_unit_weight(const unweighted_edge_list &es);

template<>
unweighted_edge_list add_unit_weight<unweighted_graph>(const unweighted_edge_list &es);

template<typename GraphType>
typename GraphType::edge_list_type add_unit_weight(const unweighted_edge_list &es) {
  typename GraphType::edge_list_type wes(es.size());
  for (size_t i = 0; i < es.size(); ++i) {
    wes[i].first = es[i].first;
    wes[i].second = typename GraphType::E{es[i].second, unit_weight<typename GraphType::W>()};
  }
  return wes;
}
}  // namespace agl
