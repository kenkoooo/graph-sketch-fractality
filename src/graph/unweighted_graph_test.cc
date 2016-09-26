#include "gtest/gtest.h"
#include "graph.h"
using namespace std;
using namespace agl;

TEST(unweighted_graph_test, instantiation) {
  G g;

  g.assign(unweighted_graph::edge_list_type{
          {0, {1}},
          {1, {2}},
          {2, {3}},
          {3, {1}},
  });

  pretty_print(g);
}

TEST(unweighted_graph_test, read_graph_tsv) {
  G g(generate_erdos_renyi(10, 2));
  pretty_print(g);
}

TEST(unweighted_graph, undirected_neighbors) {
  for (int trial = 0; trial < 10; ++trial) {
    unweighted_edge_list es = generate_erdos_renyi(100, 2);
    G g1(es);
    G g2(make_undirected(es));

    for (V v : g1.vertices()) {
      auto ns1 = range_to_vector(undirected_neighbors(g1, v));
      auto ns2 = range_to_vector<neighbor_range<unweighted_edge>>(g2.neighbors(v));
      ASSERT_TRUE(ns1 == ns2);
    }
  }
}
