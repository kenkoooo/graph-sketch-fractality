#include "generator.h"
#include <set>
#include "base/base.h"
using namespace std;

namespace agl {
unweighted_edge_list generate_path(V num_vertices) {
  unweighted_edge_list es;
  for (V v = 0; v + 1 < num_vertices; ++v) {
    es.emplace_back(v, v + 1);
  }
  return es;
}

unweighted_edge_list generate_erdos_renyi(V num_vertices, double avg_deg) {
  avg_deg = min(avg_deg, static_cast<double>(max(0, num_vertices - 1)));
  set<pair<V, V>> es;
  std::uniform_int_distribution<V> rng(0, num_vertices - 1);
  while (es.size() < num_vertices * avg_deg) {
    V u = rng(agl::random), v = rng(agl::random);
    if (u == v) continue;
    es.insert(make_pair(u, v));
  }
  return vector<pair<V, V>>(es.begin(), es.end());
};

/**
 * Generate a random scale-free network by the Barabasi-Albert (BA) model.
 * The degree distribution resulting from the BA model is scale free
 * with power-law coefficient &gamma; = 3.
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 */
unweighted_edge_list generate_ba(V final_num, V initial_num) {
  CHECK(initial_num >= 2);
  unweighted_edge_list es;
  for (int v = 0; v < initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
    }
  }

  for (int v = initial_num; v < final_num; ++v) {
    set<V> next;
    std::uniform_int_distribution<size_t> rng(0, es.size() - 1);
    while (next.size() < (size_t)initial_num) {
      size_t e = rng(agl::random);
      V u = agl::random() % 2 ? es[e].first : es[e].second;
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
    }
  }
  return es;
}

/**
 * Generate a (u, v)-flower graph by the given parameters.
 * Starting from (u + v)-vertices cycle graph, the number of vertices in the
 * resulting graph will be equal to or larger than the given number. The fractal
 * dimension D of this graph will be D = log(u + v) / log(u) (u <= v)
 * \param required_num is the minimum number of the vertices. The resulting
 * graph the number of vertices in the resulting graph will be equal to or
 * larger than this number.
 * \param u is the smaller parameter of (u, v)-flower.
 * \param v is the larger parameter of (u, v)-flower.
 */
unweighted_edge_list generate_uv_flower(V required_num, V u, V v) {
  assert(u <= v && u >= 1 && u + v >= 3);
  unweighted_edge_list es;
  es.emplace_back(u + v - 1, 0);
  for (int i = 1; i < u + v; ++i) {
    es.emplace_back(i - 1, i);
  }
  V current_vertices = u + v;
  while (current_vertices < required_num) {
    unweighted_edge_list next;
    for (auto e : es) {
      V s = e.first, t = e.second;
      for (int i = 0; i < u; ++i) {
        V left = current_vertices - 1;
        V right = current_vertices;
        if (i == 0) left = s;
        if (i == u - 1) right = t;
        next.emplace_back(left, right);
        if (right != t) current_vertices++;
      }
      for (int i = 0; i < v; ++i) {
        V left = current_vertices - 1;
        V right = current_vertices;
        if (i == 0) left = s;
        if (i == v - 1) right = t;
        next.emplace_back(left, right);
        if (right != t) current_vertices++;
      }
    }
    es.swap(next);
  }

  return es;
}

/**
 * Generate a Song-Havlin-Makse model (SHM-model) graph.
 * Starting from a star tree, the resulting graph will be a tree.
 * The fractal dimension D of this graph will be D = log(2t + 1) / log3
 * \param required_num is the minimum number of the vertices. The resulting
 * graph the number of vertices in the resulting graph will be equal to or
 * larger than this number.
 * \param initial_num is the number of vertices of the star tree of the first
 * generation.
 * \param t decides the fractal dimension of this graph. The fractal dimension D
 * of this graph will be D = log(2t + 1) / log3
 */
unweighted_edge_list generate_shm(V required_num, V initial_num, int t, double P) {
  assert(P >= 0.0 && P <= 1.0);
  std::uniform_real_distribution<> p_rng(0.0, 1.0);
  assert(t >= 2 && initial_num >= 3);
  unweighted_edge_list es;
  vector<vector<V>> adj(initial_num);
  for (int i = 1; i < initial_num; ++i) {
    es.emplace_back(0, i);
    adj[0].push_back(i);
    adj[i].push_back(0);
  }

  while ((V)adj.size() < required_num) {
    V current_num = adj.size();
    V next_num = current_num;
    for (int i = 0; i < current_num; ++i) next_num += adj[i].size() * t;

    vector<vector<V>> next(next_num);
    unweighted_edge_list next_es;
    V new_comer = current_num;
    for (V s = 0; s < current_num; ++s)
      for (int i = 0; i < (int)adj[s].size() * t; ++i) {
        next_es.emplace_back(s, new_comer);
        next[s].push_back(new_comer);
        next[new_comer].push_back(s);
        new_comer++;
      }
    for (auto e : es) {
      V s = e.first;
      V ns = -1;
      for (V n : next[s])
        if (next[n].size() == 1) {
          ns = n;
          break;
        }
      assert(ns >= 0);
      V t = e.second;
      V nt = -1;
      for (V n : next[t])
        if (next[n].size() == 1) {
          nt = n;
          break;
        }
      assert(nt >= 0);
      next_es.emplace_back(ns, nt);
      next[ns].push_back(nt);
      next[nt].push_back(ns);
    }
    if (p_rng(agl::random) <= P)
      for (auto e : es) {
        next_es.emplace_back(e);
        next[e.second].push_back(e.first);
        next[e.first].push_back(e.second);
      }
    es.swap(next_es);
    adj.swap(next);
  }
  return es;
}

unweighted_edge_list make_undirected(const unweighted_edge_list& es) {
  unweighted_edge_list out(es.size() * 2);
  for (auto i : make_irange(es.size())) {
    V u = es[i].first, v = es[i].second;
    out[i * 2 + 0] = make_pair(u, v);
    out[i * 2 + 1] = make_pair(v, u);
  }
  sort(out.begin(), out.end());
  out.erase(unique(out.begin(), out.end()), out.end());
  return out;
}

template<>
unweighted_edge_list add_random_weight<unweighted_graph>(const unweighted_edge_list &es) {
  return es;
}

template<>
unweighted_edge_list add_unit_weight<unweighted_graph>(const unweighted_edge_list &es) {
  return es;
}
}  // namespace agl
