#include "generator.h"
#include "gtest/gtest.h"
#include <vector>
#include <random>
#include <queue>

using namespace agl;
using namespace std;

namespace {
template<typename GraphType>
bool is_connected(const GraphType &g) {
  if (g.num_vertices() == 0) return true;

  std::queue<V> que;
  std::vector<bool> vis(g.num_vertices());
  que.push(0);
  vis[0] = true;
  while (!que.empty()) {
    V v = que.front();
    que.pop();
    for (V tv : g.neighbors(v)) {
      if (vis[tv]) continue;
      que.push(tv);
      vis[tv] = true;
    }
  }
  return std::find(vis.begin(), vis.end(), false) == vis.end();
}
}  // namespace

TEST(gen_ba, random_num_vertices) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = agl::random(1000) + 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_ba, corner_case) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 2;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_ba, check_double_edge) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 2;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    set<pair<V, V>> s;
    for (const auto &e : es) s.insert(e);
    ASSERT_EQ(s.size(), es.size());

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_uv_flower, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V u = agl::random(3) + 1;
    V v = u + agl::random(3) + 1;
    V w = u + v;
    V req = w;
    int n = agl::random(3);
    for (int i = 0; i < n; ++i) {
      req *= w;
    }
    cerr << req << " " << u << " " << v << endl;
    auto es = generate_uv_flower(req, u, v);

    // Number of edges
    size_t expected_edge_num = w;
    size_t expected_node_num = w;
    size_t max_deg = 2;
    while (expected_edge_num < es.size()) {
      expected_edge_num *= w;
      expected_node_num = w * expected_node_num - w;
      max_deg *= 2;
    }
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(make_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
    ASSERT_EQ(g.num_vertices(), expected_node_num);

    vector<V> deg_check(expected_node_num, 0);
    V current = w;
    for (int i = 0; i < current; ++i) deg_check[i] = 2;
    while (current < expected_node_num) {
      current = current * w - w;
      for (int i = 0; i < current; ++i) {
        if (deg_check[i] == 0)
          deg_check[i] = 2;
        else
          deg_check[i] *= 2;
      }
    }

    // Check degree
    for (int i = 0; i < g.num_vertices(); ++i) {
      ASSERT_EQ(g.degree(i), deg_check[i]);
    }
  }
}

TEST(gen_shm, small_case) {
  V initial_num = 5;
  V required_num = initial_num;
  int t = 2;
  int generation = 3;
  size_t max_deg = initial_num - 1;
  for (int i = 1; i < generation; ++i) {
    required_num = (2 * t + 1) * required_num - 2 * t;
    max_deg *= t;
  }
  auto es = generate_shm(required_num, initial_num, t);

  // Number of edges
  ASSERT_EQ(es.size(), required_num - 1);

  G g(make_undirected(es));
  pretty_print(g);
  ASSERT_TRUE(is_connected(g));
  ASSERT_EQ(g.num_vertices(), required_num);
  ASSERT_EQ(g.degree(0), max_deg);

  // Check degree
  for (V v : g.vertices()) {
    V deg = g.degree(v);
    while (deg % t == 0) deg /= t;
    ASSERT_TRUE(deg == initial_num - 1 || deg == 1 || deg == 2);
  }
}

TEST(gen_shm, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V initial_num = agl::random(20) + 3;
    V required_num = initial_num;
    int t = agl::random(20) + 2;
    int generation = agl::random(5) + 1;
    size_t max_deg = initial_num - 1;
    for (int i = 1; i < generation; ++i) {
      required_num = (2 * t + 1) * required_num - 2 * t;
      max_deg *= t;
    }
    auto es = generate_shm(required_num, initial_num, t);

    // Number of edges
    ASSERT_EQ(es.size(), required_num - 1);

    G g(make_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
    ASSERT_EQ(g.num_vertices(), required_num);
    ASSERT_EQ(g.degree(0), max_deg);
    // Check degree
    for (V v : g.vertices()) {
      if (v == 0) continue;
      V deg = g.degree(v);
      while (deg % t == 0) deg /= t;
      ASSERT_TRUE(deg == initial_num - 1 || deg == 1 || deg == 2);
    }
  }
}

TEST(gen_shm, small_world) {
  for (int trial = 0; trial < 10; ++trial) {
    V initial_num = 5;
    V required_num = initial_num;
    int t = 2;
    int generation = agl::random(5) + 1;
    size_t max_deg = initial_num - 1;
    size_t expected_edge_num = initial_num - 1;
    for (int i = 1; i < generation; ++i) {
      required_num = required_num + 2 * t * expected_edge_num;
      expected_edge_num = expected_edge_num * (2 * t + 2);
      max_deg *= t + 1;
    }
    auto es = generate_shm(required_num, initial_num, t, 1.0);

    // Number of edges
    ASSERT_EQ(es.size(), expected_edge_num);
    G g(make_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
    ASSERT_EQ(g.num_vertices(), required_num);
    ASSERT_EQ(g.degree(0), max_deg);
  }
}
