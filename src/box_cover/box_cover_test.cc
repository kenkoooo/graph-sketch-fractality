#include "box_cover.h"
#include "gtest/gtest.h"

using namespace agl;
using namespace agl::box_cover_internal;
using namespace std;

pair<vector<V>, vector<V>> inv_and_rank(const G &g) {
  vector<V> rank(g.num_vertices());
  vector<V> inv(g.num_vertices());
  for (V i = 0; i < g.num_vertices(); ++i) inv[i] = i;
  shuffle(inv.begin(), inv.end(), agl::random);
  for (int i = 0; i < g.num_vertices(); ++i) rank[inv[i]] = i;
  return {inv, rank};
}

TEST(box_cover, memb) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = 3;
    vector<V> memb = box_cover_memb(g, radius);

    coverage_manager cm(g, radius, 1.0);
    for (const V &v : memb) cm.add(g, v);

    ASSERT_EQ(cm.get_current_coverage(), 1.0);
  }
}

TEST(box_cover, burning) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(200);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = 1;
    vector<V> burning = box_cover_burning(g, radius);

    coverage_manager cm(g, radius, 1.0);
    for (const V &v : burning) cm.add(g, v);

    ASSERT_EQ(cm.get_current_coverage(), 1.0);
  }
}

TEST(box_cover, build_sketch_check) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = agl::random(3) + 1;
    const int k = 128;

    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;

    coverage_manager cm(g, radius, 1.0);
    for (int cover_trial = 0; cover_trial < 10; ++cover_trial) {
      vector<bool> covered(g.num_vertices());
      for (int i = 0; i < g.num_vertices(); ++i) covered[i] = cm.v_covered(i);

      // build sketches
      vector<vector<V>> naive_x =
          naive_build_sketch(g, radius, k, rank, inv, covered);
      vector<vector<V>> x = build_sketch(g, radius, k, rank, inv, cm);

      // check
      for (V v = 0; v < g.num_vertices(); v++) ASSERT_EQ(naive_x[v], x[v]) << v;
    }
  }
}

TEST(box_cover, greedy_small) {
  for (int trial = 0; trial < 1000; ++trial) {
    const W radius = agl::random(3) + 1;
    V M = 3;
    V N = M + agl::random(50);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    const int k = agl::random(20) + 5;
    if (g.num_vertices() <= k) {
      trial--;
      continue;
    }
    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;

    coverage_manager cm(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);

    vector<V> centers2;
    {
      vector<bool> centered(g.num_vertices(), false);
      naive_select_greedily(g, X, centers2, centered, k);
    }

    vector<V> centers1;
    {
      coverage_manager cm(g, radius, 1.0);
      select_greedily(g, X, centers1, k, cm);
    }
    ASSERT_EQ(centers1, centers2);
  }
}

TEST(box_cover, greedy_big) {
  for (int trial = 0; trial < 30; ++trial) {
    const W radius = agl::random(4) + 1;
    V M = 3;
    V N = M + agl::random(2000) + 2000;
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    const int k = 1024;
    if (g.num_vertices() < k) {
      trial--;
      continue;
    }
    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;

    coverage_manager cm(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);
    vector<V> centers1;
    {
      coverage_manager cm(g, radius, 1.0);
      select_greedily(g, X, centers1, k, cm);
    }
    vector<V> centers2;
    {
      vector<bool> centered(g.num_vertices(), false);
      naive_select_greedily(g, X, centers2, centered, k);
    }
    ASSERT_EQ(centers1, centers2);
    cerr << "Test case " << (trial + 1) << " has been done." << endl;
  }
}

TEST(box_cover, greedy_huge) {
  const int k = 1024;
  W radius = agl::random(4) + 1;
  V M = 3;
  V N = M + agl::random(10000) + 10000;
  while (N <= k) N = M + agl::random(10000) + 10000;
  auto es = generate_ba(N, M);
  G g(make_undirected(es));
  auto inv_rank = inv_and_rank(g);
  auto &inv = inv_rank.first;
  auto &rank = inv_rank.second;

  coverage_manager cm(g, radius, 1.0);
  vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);
  pretty_print(g);
  vector<V> centers1;
  {
    coverage_manager cm(g, radius, 1.0);
    select_greedily(g, X, centers1, k, cm);
  }
  vector<V> centers2;
  {
    vector<bool> centered(g.num_vertices(), false);
    naive_select_greedily(g, X, centers2, centered, k);
  }

  ASSERT_EQ(centers1, centers2);
}

TEST(box_cover, coverage_management) {
  for (int trial = 0; trial < 100; ++trial) {
    const int k = 1024;
    W radius = agl::random(2) + 1;
    V M = 3;
    V N = M + agl::random(1000) + 1000;
    while (N <= k) N = M + agl::random(1000) + 1000;
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;

    coverage_manager cmtmp(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cmtmp);
    pretty_print(g);
    vector<V> centers;
    coverage_manager cm(g, radius, 1.0);
    select_greedily(g, X, centers, k, cm);
    double tester = naive_coverage(g, centers, radius);
    double hey = cm.get_current_coverage();
    ASSERT_EQ(tester, hey);
  }
}

TEST(box_cover, coverage_break) {
  for (int trial = 0; trial < 100; ++trial) {
    const int k = 1024;
    W radius = agl::random(2) + 1;
    V M = 3;
    V N = M + agl::random(1000) + 1000;
    while (N <= k) N = M + agl::random(1000) + 1000;
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;

    coverage_manager cmtmp(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cmtmp);
    pretty_print(g);
    vector<bool> centered(g.num_vertices(), false);

    double goal = (double)(agl::random(5) + 95) / 100;
    vector<V> centers = box_cover_sketch(g, radius, k, 100, goal);

    coverage_manager cm(g, radius, goal);
    for (V c : centers) cm.add(g, c);
    ASSERT_TRUE(cm.get_current_coverage() >= goal);
  }
}

TEST(box_cover, find_solution_flower) {
  vector<pair<V, V>> uvs = {{1, 2}, {2, 2}};
  for (auto p : uvs) {
    V u = p.first, v = p.second;
    V req = 10000;
    auto es = generate_uv_flower(req, u, v);
    G g(make_undirected(es));
    auto pairs = find_analytical_solution("flower", u, v, g);
    for (auto p : pairs) {
      cerr << p.first << " " << p.second << endl;
    }
  }
}

TEST(box_cover, find_solution_shm) {
  vector<pair<V, int>> nts = {{5, 2}};
  for (auto p : nts) {
    V n = p.first, t = p.second;
    V req = 50;
    auto es = generate_shm(req, n, t);
    G g(make_undirected(es));
    auto pairs = find_analytical_solution("shm", n, t, g);
    for (auto p : pairs) {
      cerr << p.first << " " << p.second << endl;
    }
  }
}

TEST(box_cover, lazy_greedily) {
  V u = 2, v = 2;
  auto es = generate_uv_flower(10000, u, v);
  G g(make_undirected(es));
  vector<pair<W, V>> pairs = find_analytical_solution("flower", u, v, g);

  const int k = 1024;
  for (auto p : pairs) {
    W rad = p.first / 2;

    coverage_manager cm(g, rad, 1.0);
    bool tmp = true;
    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;
    auto X = build_sketch(g, rad, k, rank, inv, cm, tmp, SIZE_MAX);

    vector<V> centers;
    select_lazy_greedily(g, X, rank, inv, centers, cm);

    vector<V> memb = box_cover_memb(g, rad);

    size_t memb_size = memb.size();
    size_t sketch_size = centers.size();
    ASSERT_TRUE(memb_size < sketch_size * 2 && sketch_size < memb_size * 2);
  }
}

TEST(box_cover, coloring) {
  V u = 2, v = 2;
  auto es = generate_uv_flower(1000, u, v);
  G g(make_undirected(es));
  pretty_print(g);

  vector<pair<W, V>> pairs = find_analytical_solution("flower", u, v, g);
  map<W, size_t> m;
  for (const auto &p : pairs) m[p.first] = p.second;

  W largest = pairs.rbegin()->first;
  vector<pair<W, size_t>> coloring = box_cover_coloring(g, largest);
  for (auto c : coloring)
    if (m[c.first] > 0)
      ASSERT_TRUE(c.second < m[c.first] * 10 && m[c.first] < c.second * 10);
}

TEST(box_cover, covered_check) {
  for (int trial = 0; trial < 10; ++trial) {
    V N = agl::random(1000) + 1000;
    W radius = agl::random(10) + 10;
    int k = 128;
    auto es = generate_uv_flower(N, 2, 2);
    G g(make_undirected(es));
    coverage_manager cm(g, radius, 1.0);
    auto inv_rank = inv_and_rank(g);
    auto &inv = inv_rank.first;
    auto &rank = inv_rank.second;

    cm.add(g, 0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);
    ASSERT_TRUE(X[0].empty());
  }
}

TEST(box_cover, cbb) {
  for (int trial = 0; trial < 10; ++trial) {
    auto es = generate_uv_flower(agl::random(1000) + 1000, 2, 2);
    G g(make_undirected(es));
    W d = 2 + agl::random(10);

    auto cbb = box_cover_cbb(g, d);
    auto memb = box_cover_memb(g, d / 2);
    ASSERT_TRUE(cbb.size() * 10 > memb.size() && memb.size() * 10 > cbb.size());
  }
}
