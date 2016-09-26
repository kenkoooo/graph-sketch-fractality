#include "agl.h"
#include <vector>
#include <set>
#include <map>
#include <queue>

namespace agl {
namespace box_cover_internal {
class coverage_manager {
 private:
  V num_v;
  std::vector<W> dist;
  V cnt;
  W radius;
  double least_coverage;
  coverage_manager() {}

 public:
  /**
   * Coverage Manager
   * \param g is the graph to cover
   * \param r is the radius of each box
   * \param c is the least coverage
   */
  coverage_manager(const G &g, const W r, const double c)
      : cnt(0), radius(r), least_coverage(c) {
    assert(least_coverage <= 1.0);
    num_v = g.num_vertices();
    dist.assign(num_v, num_v);
  }

  void add(const G &g, const V new_v) {
    std::queue<V> que;
    if (dist[new_v] == num_v) cnt++;
    dist[new_v] = 0;
    que.push(new_v);

    for (W d = 0; d < radius; ++d) {
      size_t s = que.size();
      for (size_t t = 0; t < s; t++) {
        V a = que.front();
        que.pop();
        for (V neighbor : g.neighbors(a)) {
          if (dist[neighbor] <= d + 1) continue;
          if (dist[neighbor] == num_v) cnt++;
          dist[neighbor] = d + 1;
          que.push(neighbor);
        }
      }
    }
  }
  double get_current_coverage() { return (double)cnt / num_v; }
  bool is_covered() { return get_current_coverage() >= least_coverage; }
  bool v_covered(V v) const { return dist[v] <= radius; }
  bool is_center(V v) const { return dist[v] == 0; }
};

std::vector<std::pair<W, V>> find_analytical_solution(const std::string &type, V u, V v, const G &g);
void select_greedily(const G &g, const std::vector<std::vector<V>> &X, std::vector<V> &centers, const int k, coverage_manager &cm);
void select_lazy_greedily(const G &g, const std::vector<std::vector<V>> &X, const std::vector<V> &rank, const std::vector<V> &inv, std::vector<V> &centers, coverage_manager &cm);
std::vector<std::vector<V>> build_sketch(const G &g, const W radius, const int k, const std::vector<V> &rank, const std::vector<V> &inv, const coverage_manager &cm);
std::vector<std::vector<V>> build_sketch(const G &g, const W radius, const int k, const std::vector<V> &rank, const std::vector<V> &inv, const coverage_manager &cm, bool &use_memb, size_t index_size_limit);

//
// Naive Functions for Tests
//
double naive_coverage(const G &g, const std::vector<V> &s, W rad);
std::vector<std::vector<V>> naive_build_sketch(const G &g, const W radius, const int k, const std::vector<V> &rank, const std::vector<V> &inv, const std::vector<bool> &is_covered);
void naive_select_greedily(const G &g, const std::vector<std::vector<V>> &X, std::vector<V> &centers, std::vector<bool> &centered, const int k);

}  // namespace box_cover_internal

//
// Radius-based Methods:
//   Returns the set S of selected center nodes.
//   For any vertex v, there is s \in S s.t. d(v, s) <= radius.
//

//! Song et al. 2007 (Section 3.2)
std::vector<V> box_cover_memb(const G &g, W radius);

//! Schneider et al. 2012
std::vector<V> box_cover_burning(const G &g, W radius);

//! Akiba et al. 2016
std::vector<V> box_cover_sketch(const G &g, W radius, const int k, const int pass, double least_coverage = 1.0, double alpha = 1.0);
std::vector<V> box_cover_sketch(const G &g, W radius, const int k, const int pass, box_cover_internal::coverage_manager &cm, double alpha = 1.0);

//
// Diameter-based Methods:
//   returns the sets of vertices with the limited diameter.
//   Any vertex is covered by a set.
//

//! Song et al. 2007 (Section 3.1)
std::vector<V> box_cover_cbb(const G &g, W diameter);

//! Song et al. 2007 (Section 2)
std::vector<std::pair<W, size_t>> box_cover_coloring(const G &g, W diameter);
}  // namespace agl
