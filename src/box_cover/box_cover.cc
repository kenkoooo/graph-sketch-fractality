#include "box_cover.h"
#include <queue>
using namespace std;

namespace agl{
namespace box_cover_internal {

double naive_coverage(const G &g, const vector<V> &s, W rad, vector<bool> &is_covered) {
  vector<W> dist(g.num_vertices(), g.num_vertices());
  for (const V &start : s) {
    queue<pair<V, W>> que;
    que.push({start, 0});
    is_covered[start] = true;
    while (!que.empty()) {
      V v = que.front().first;
      W d = que.front().second;
      que.pop();
      dist[v] = d;
      if (d == rad) continue;
      for (const V &n : g.neighbors(v)) {
        if (dist[n] < d + 1) continue;
        is_covered[n] = true;
        que.push({n, d + 1});
      }
    }
  }
  V cnt = 0;
  for (const bool &c : is_covered)
    if (c) cnt++;
  return (double)cnt / g.num_vertices();
}

double naive_coverage(const G &g, const vector<V> &s, W rad) {
  vector<bool> is_covered_dummy(g.num_vertices(), false);
  return naive_coverage(g, s, rad, is_covered_dummy);
}

vector<V> merge_and_purify(set<V> &parent, const vector<V> &sorted_vec, const int k) {
  vector<V> delta;
  for (const V &p_rank : sorted_vec) {
    if (parent.size() == (size_t)k && p_rank > *parent.rbegin()) break;
    size_t prev = parent.size();
    parent.insert(p_rank);
    if (parent.size() > prev) delta.push_back(p_rank);
    while (parent.size() > (size_t)k) parent.erase(*parent.rbegin());
  }
  return delta;
}

/**
 * Naive method to generate bottom-k min-hash sketch.
 * \param g is graph to cover
 * \param radius is radius of each box
 * \param k size of each bottom-k min-hash sketch
 * \param rank Ranks of vertices.
 * \param inv Inverted index of rank.
 * \param is_covered Is vertex v is covered?
 */
vector<vector<V>> naive_build_sketch(const G &g, const W radius, const int k, const vector<V> &rank, const vector<V> &inv, const vector<bool> &is_covered) {
  V num_v = g.num_vertices();
  vector<vector<V>> naive_X(num_v);
  for (V i = 0; i < num_v; ++i) {
    set<V> tmp;
    vector<bool> vis(num_v, false);
    queue<pair<V, W>> que;
    que.push(make_pair(i, 0));
    vis[i] = true;
    tmp.insert(rank[i]);

    while (!que.empty()) {
      V v = que.front().first;
      W dist = que.front().second;
      que.pop();
      if (dist == radius) break;
      for (V u : g.neighbors(v)) {
        if (vis[u]) continue;
        que.push(make_pair(u, dist + 1));
        vis[u] = true;
        tmp.insert(rank[u]);
      }
    }
    int cnt = 0;
    for (const V &p : tmp) {
      if (is_covered[inv[p]]) continue;
      naive_X[i].push_back(p);
      cnt++;
      if (cnt == k) break;
    }
  }

  return naive_X;
}

/**
 * Generate bottom-k min-hash sketch.
 * \param g the graph to cover
 * \param radius radius of each box
 * \param k maximum size of each sketch
 * \param rank Ranks of vertices.
 * \param inv Inverted index of rank.
 * \param cm coverage manager
 */
vector<vector<V>> build_sketch(const G &g, const W radius, const int k, const vector<V> &rank, const vector<V> &inv, const coverage_manager &cm) {
  bool t = false;
  return build_sketch(g, radius, k, rank, inv, cm, t, g.num_vertices() * k);
}

/**
 * Generate bottom-k min-hash sketch.
 * \param g the graph to cover
 * \param radius radius of each box
 * \param k maximum size of each sketch
 * \param rank Ranks of vertices.
 * \param inv Inverted index of rank.
 * \param use_memb Want to use MEMB if the required memory is enough small?
 * \param index_size_limit limit of memory size which MEMB can be used.
 */
vector<vector<V>> build_sketch(const G &g, const W radius, const int k, const vector<V> &rank, const vector<V> &inv, const coverage_manager &cm, bool &use_memb, size_t index_size_limit) {
  V num_v = g.num_vertices();
  vector<set<V>> X(num_v);
  vector<vector<V>> previous_added(num_v);

  //
  // Build-Sketches O((n+m)*rad)
  //
  size_t using_k = use_memb ? X[0].max_size() - 1 : k;
  size_t total_size = 0;

  for (V i = 0; i < num_v; ++i) {
    if (cm.v_covered(i)) continue;
    X[i].insert(rank[i]);
    previous_added[i].push_back(i);
    total_size++;
  }

  for (W d = 0; d < radius; ++d) {
    for (const V &v : inv) {
      vector<V> next;
      for (const V &a : previous_added[v]) {
        for (const V &neighbor : g.neighbors(a)) {
          // Merge & Purify
          V rv = rank[v];
          set<V> &xn = X[neighbor];

          if (xn.size() >= using_k && *(xn.rbegin()) <= rv) continue;

          auto inserted = xn.insert(rv);
          while (xn.size() > using_k) {
            auto it = xn.end();
            it--;
            xn.erase(it);
          }
          if (inserted.second && *xn.rbegin() >= rv) {
            next.push_back(neighbor);
            total_size++;
          }
        }

        if (use_memb && total_size >= index_size_limit) {
          using_k = k;
          use_memb = false;
          for (V i = 0; i < num_v; ++i)
            while (X[i].size() > using_k) {
              auto it = X[i].end();
              it--;
              X[i].erase(it);
            }
        }
      }
      previous_added[v].swap(next);
    }
  }

  vector<vector<V>> ret;
  for (V v = 0; v < num_v; ++v) {
    while (X[v].size() > using_k) X[v].erase(*(X[v].rbegin()));
    vector<V> sketch(X[v].begin(), X[v].end());
    ret.push_back(sketch);
  }
  return ret;
}

/**
 * Select center nodes of boxes by MEMB
 * \param g is graph to cover
 * \param X bottom-k min-hash sketch.
 * \param rank Ranks of vertices.
 * \param inv Inverted index of rank.
 * \param centers center nodes of boxes
 * \param cm coverage manager
 */
void select_lazy_greedily(const G &g, const vector<vector<V>> &X, const vector<V> &rank, const vector<V> &inv, vector<V> &centers, coverage_manager &cm) {
  vector<bool> rank_covered(X.size(), false);
  priority_queue<pair<V, V>> que;
  vector<V> box_size(X.size());
  for (size_t box = 0; box < X.size(); ++box) {
    if (X[box].size() == 0 || cm.is_center(box)) continue;
    box_size[box] = X[box].size();
    que.push({box_size[box], rank[box]});
  }

  vector<vector<V>> inverted(X.size());
  for (size_t box = 0; box < X.size(); ++box)
    for (const V &rank_b : X[box]) inverted[rank_b].push_back(box);

  while (!que.empty() && !cm.is_covered()) {
    V s = que.top().first;
    V rank_box = que.top().second;
    V v = inv[rank_box];
    que.pop();
    if (cm.is_center(v)) continue;
    if (box_size[v] != s) {
      que.push({box_size[v], rank[v]});
      continue;
    }

    centers.push_back(v);
    cm.add(g, v);
    if (cm.is_covered()) break;

    for (const auto &rank_v : X[v]) {
      if (rank_covered[rank_v]) continue;
      rank_covered[rank_v] = true;
      for (const auto &box : inverted[rank_v]) box_size[box]--;
    }
  }
}

/**
 * Select center nodes of boxes by using min-hash sketch
 * \param g is graph to cover
 * \param X bottom-k min-hash sketch.
 * \param centers center nodes of boxes
 * \param k maximum size of each sketch
 * \param cm coverage manager
 */
void select_greedily(const G &g, const vector<vector<V>> &X, vector<V> &centers, const int k, coverage_manager &cm) {
  assert(g.num_vertices() > k);
  //
  // Variables
  //
  V num_v = g.num_vertices();
  set<V> Xs;
  priority_queue<pair<V, V>, vector<pair<V, V>>, greater<pair<V, V>>> que[2];
  vector<multimap<V, V>> T(k + 2);
  vector<V> k1(num_v);
  vector<V> k2(num_v);
  vector<bool> is_type1(num_v, false);
  vector<vector<V>> I(num_v);
  vector<bool> removed(num_v, false);
  vector<bool> covered_rank(num_v, false);

  auto last_element = [&](V box) -> V { return X[box][k1[box] - 1]; };
  auto insert_as_type = [&](V box, int target_type) {
    if (target_type == 0) {
      is_type1[box] = false;
      T[k2[box] + 1].insert({last_element(box), box});
      que[0].push({last_element(box), box});
    } else {
      is_type1[box] = true;
      T[k2[box]].insert({last_element(box), box});
      que[1].push({k2[box], box});
    }
  };
  auto remove_covered_ranks = [&](V box) {
    while (covered_rank[last_element(box)]) {
      k1[box]--;
      if (k1[box] == 0) {
        removed[box] = true;
        break;
      }
    }
  };
  auto remove_multimap_pair = [&](V from, const pair<V, V> &p) {
    auto it_p = T[from].equal_range(p.first);
    for (auto it = it_p.first; it != T[from].end(); it++) {
      if ((*it).second == p.second) {
        T[from].erase(it);
        break;
      }
      if (it == it_p.second) break;
    }
  };

  //
  // Initialization
  //
  for (V p = 0; p < num_v; ++p) {
    if (cm.is_center(p) || X[p].empty()) continue;
    k1[p] = X[p].size();
    k2[p] = k - k1[p];
    if (X[p].size() == (size_t)k) {
      que[0].push({last_element(p), p});
      is_type1[p] = false;
      T[k2[p] + 1].insert({last_element(p), p});
    } else {
      que[0].push({num_v, p});
      que[1].push({k2[p], p});
      is_type1[p] = true;
      T[k2[p]].insert({last_element(p), p});
    }
    for (V ri : X[p]) {
      I[ri].push_back(p);
    }
  }

  //
  // Main loop
  //
  while ((Xs.size() == (size_t)k ? *Xs.rbegin() : num_v) > k - 1) {
    // Selection
    V select = -1;
    V argmin = num_v + 1;
    for (int q = 0; q < 2; q++) {
      while (!que[q].empty()) {
        // Remove unnecessary or changed elements from queue
        V top_key = que[q].top().first;
        V top_v = que[q].top().second;

        if (top_key == num_v && Xs.size() < (size_t)k && !cm.is_center(top_v)) {
          break;  // Fit to Naive Method
        } else if (removed[top_v] || q != is_type1[top_v]) {
          que[q].pop();
          continue;
        } else if (q == 1 && k2[top_v] != top_key) {
          que[q].pop();
          continue;
        }
        break;
      }
      if (que[q].empty()) continue;
      V v = que[q].top().second;

      set<V> tmp(Xs);
      merge_and_purify(tmp, X[v], k);
      V ec_tmp = tmp.size() == (size_t)k ? *tmp.rbegin() : num_v;
      if (argmin > ec_tmp || (argmin == ec_tmp && v < select)) {
        argmin = ec_tmp;
        select = v;
      }
    }
    if (select < 0) break;
    assert(!cm.is_center(select));

    //
    // Merge and Remove
    //
    vector<V> delta = merge_and_purify(Xs, X[select], k);
    centers.push_back(select);
    removed[select] = true;
    remove_multimap_pair(k2[select], {last_element(select), select});
    cm.add(g, select);
    if (cm.is_covered()) return;

    //
    // Update about covered elements
    //
    for (const V &rank_i : delta) {
      covered_rank[rank_i] = true;
      for (const V &box : I[rank_i]) {
        if (k2[box] + 1 >= k) removed[box] = true;
        if (removed[box]) continue;
        pair<V, V> box_pair = {last_element(box), box};

        // When the subbox is Type0 and its last element is covered,
        // it has to be Type1

        if (is_type1[box]) {
          remove_multimap_pair(k2[box], box_pair);
          T[k2[box] + 1].insert(box_pair);
          que[1].push({k2[box] + 1, box});
        } else if (rank_i == last_element(box)) {
          remove_multimap_pair(k2[box] + 1, box_pair);
          remove_covered_ranks(box);
          if (removed[box]) continue;
          is_type1[box] = true;
          if (k2[box] + 1 <= k) T[k2[box] + 1].insert({last_element(box), box});
          que[1].push({k2[box] + 1, box});
        } else {
          remove_multimap_pair(k2[box] + 1, box_pair);
          T[k2[box] + 2].insert(box_pair);
        }
        k2[box]++;
      }
    }

    //
    // Update
    //
    int j = 1;
    for (auto it_Xs = Xs.begin(); it_Xs != Xs.end(); ++it_Xs, ++j) {
      if (T[j].empty()) continue;

      vector<pair<V, V>> removing;  // removed elements from T after iteration
      for (auto it = T[j].lower_bound(*it_Xs); it != T[j].end(); ++it) {
        V box = (*it).second;

        auto it_j = it_Xs;
        // always last_blue[box] >= jth_rank
        removing.push_back({last_element(box), box});
        if (removed[box]) continue;
        if (!is_type1[box]) {  // Type0->
          assert(k2[box] + 1 == j);
          if (covered_rank[last_element(box)]) {
            k2[box]--;
            it_j--;
          }
          k1[box]--;
          k2[box]++;
          if (k1[box] == 0 || k2[box] >= k) {
            removed[box] = true;
            continue;
          }
        }

        remove_covered_ranks(box);
        if (removed[box]) continue;
        if (last_element(box) > *it_j) {  // Type0->Type0
          insert_as_type(box, 0);
        } else {  // Type0->Type1
          insert_as_type(box, 1);
        }
      }
      for (auto rm_pair : removing) remove_multimap_pair(j, rm_pair);
    }
  }
}

/**
 * naive method to select center nodes of boxes by using min-hash sketch
 * \param g is graph to cover
 * \param X bottom-k min-hash sketch.
 * \param centers center nodes of boxes
 * \param centered Is node v center?
 * \param k k
 */
void naive_select_greedily(const G &g, const vector<vector<V>> &X, vector<V> &centers, vector<bool> &centered, const int k) {
  V num_v = g.num_vertices();
  set<V> Xs;
  auto estimated_cardinality = [&](const set<V> &subset) -> double {
    if (subset.size() < (size_t)k) {
      return (k - 1);
    }
    if (subset.size() == (size_t)k) {
      V kth = *subset.rbegin();
      return (double)(k - 1) / kth * g.num_vertices();
    }
    assert(false);
    return (k - 1);
  };

  while (estimated_cardinality(Xs) < max(num_v, k - 1)) {
    V selected_v = -1;
    double argmax = 0.0;
    for (V v = 0; v < num_v; v++) {
      if (centered[v]) continue;
      set<V> tmp(Xs);
      merge_and_purify(tmp, X[v], k);
      double ec_tmp = estimated_cardinality(tmp);
      if (argmax < ec_tmp) {
        argmax = ec_tmp;
        selected_v = v;
      }
    }

    if (selected_v < 0) {
      break;
    }

    centers.push_back(selected_v);
    centered[selected_v] = true;

    merge_and_purify(Xs, X[selected_v], k);
  }
}

/**
 * Calculate analytical size of boxes for two fractal models, (u, v)-flower and
 * SHM-model.
 * \param type is model type, which has to be "flower" or "shm"
 * \param u parameter u of (u, v)-flower, or initial number of nodes of
 * SHM-model
 * \param v parameter v of (u, v)-flower, or multiply parameter of SHM-model
 * \param g graph to calculate
 */
vector<pair<W, V>> find_analytical_solution(const string &type, V u, V v, const G &g) {
  vector<W> diameters;
  vector<V> nodes;
  vector<V> edges;
  if (type == "flower") {
    if (u > v) swap(u, v);
    V w = u + v;
    V M = w;
    V N = u + v;
    W d = (u + v) / 2;
    edges.push_back(M);
    nodes.push_back(N);
    diameters.push_back(d);
    while (N < g.num_vertices()) {
      N = w * N - w;
      d = u * d + (v - u);
      M *= w;
      edges.push_back(M);
      nodes.push_back(N);
      diameters.push_back(d);
    }
    assert(N == g.num_vertices() && M == (V)g.num_edges() / 2);

    int n = edges.size();
    vector<pair<W, V>> ret;
    for (int m = 1; m <= n; ++m) {
      V Nb = (w - 2) * pow(w, n - m) + w;
      assert(Nb % (w - 1) == 0);
      Nb /= (w - 1);
      ret.push_back({diameters[m - 1], Nb});
    }
    return ret;
  } else if (type == "shm") {
    W d = 2;
    V N = u;
    int t = v;
    V M = u - 1;
    nodes.push_back(1);
    nodes.push_back(N);
    edges.push_back(M);
    diameters.push_back(d);
    while (N < g.num_vertices()) {
      N = N + 2 * t * M;
      M = N - 1;
      d = 3 * d + 2;
      edges.push_back(M);
      nodes.push_back(N);
      diameters.push_back(d);
    }
    assert(N == g.num_vertices() && M == (V)g.num_edges() / 2);
    int n = edges.size();
    vector<pair<W, V>> ret;
    for (int m = 1; m <= n; ++m) {
      ret.push_back({diameters[m - 1], nodes[n - m]});
    }
    return ret;

  } else {
    cerr << "Unknown graph type: " << type << endl;
  }
  return {};
}
}  // namespace box_cover_internal

using namespace agl::box_cover_internal;

/**
 * Box-Covering algorithm, named maximum excluded mass burning (MEMB).
 * \param g is graph to cover
 * \param radius is radius of each box
 */
vector<V> box_cover_memb(const G &g, W radius) {
  V num_v = g.num_vertices();
  vector<vector<pair<V, W>>> node_lists;
  map<size_t, set<V>> excluded_mass_map;
  {
    vector<pair<size_t, V>> center_candidates;
    for (V pv = 0; pv < num_v; ++pv) {
      queue<pair<V, W>> que;
      vector<bool> vis(num_v, false);
      vector<pair<V, W>> nodes;

      que.push(make_pair(pv, 0));
      nodes.push_back(make_pair(pv, 0));
      vis[pv] = true;

      while (!que.empty()) {
        V v = que.front().first;
        W dist = que.front().second;
        que.pop();
        if (dist >= radius) continue;
        for (const V &u : g.neighbors(v)) {
          if (vis[u]) continue;
          que.push(make_pair(u, dist + 1));
          nodes.push_back(make_pair(u, dist + 1));
          vis[u] = true;
        }
      }

      node_lists.push_back(nodes);
      center_candidates.push_back(make_pair(nodes.size(), pv));
    }
    for (const pair<size_t, V> &p : center_candidates) {
      excluded_mass_map[p.first].insert(p.second);
    }
  }

  set<V> covered_nodes;
  set<V> center_nodes;
  vector<W> central_distance(num_v, num_v);
  while (covered_nodes.size() < (size_t)num_v) {
    V center_node_found;
    while (true) {
      V node, maximum_key;
      while (true) {
        maximum_key = excluded_mass_map.rbegin()->first;
        set<V> &nodes = excluded_mass_map.rbegin()->second;
        auto it = nodes.begin();
        advance(it, agl::random(nodes.size()));
        node = *it;
        if (center_nodes.find(node) != center_nodes.end()) {
          nodes.erase(it);
          if (nodes.empty()) excluded_mass_map.erase(maximum_key);
        } else {
          break;
        }
      }

      V mass = 0;
      for (const pair<V, W> &p : node_lists[node]) {
        if (covered_nodes.find(p.first) == covered_nodes.end()) {
          mass++;
        }
      }
      excluded_mass_map[maximum_key].erase(node);
      if (excluded_mass_map[maximum_key].empty())
        excluded_mass_map.erase(maximum_key);
      if (mass == maximum_key) {
        center_node_found = node;
        break;
      } else {
        excluded_mass_map[mass].insert(node);
      }
    }
    center_nodes.insert(center_node_found);
    for (const pair<V, W> &p : node_lists[center_node_found]) {
      V i = p.first;
      W d = p.second;
      covered_nodes.insert(i);
      central_distance[i] = max(central_distance[i], d);
    }
  }

  vector<V> ret(center_nodes.begin(), center_nodes.end());
  return ret;
}

/**
 * Box-Covering algorithm, introduced by Schneider et al. in 2012.
 * The box-covering for tree networks could be performed in O(N^3) while for
 * regular networks it requires O(2^N).
 * \param g is graph to cover
 * \param radius is radius of each box
 */
vector<V> box_cover_burning(const G &g, W radius) {
  if (radius == 0) {
    vector<V> ret;
    for (V v = 0; v < g.num_vertices(); v++) ret.push_back(v);
    return ret;
  }
  auto burning_splitted = [](const G &g, W radius, set<V> &solution,
                             vector<vector<V>> &boxes) {
    V num_v = g.num_vertices();
    V prev_size = -1;
    while (prev_size < (V)solution.size()) {
      prev_size = solution.size();
      //
      // Remove unnecessary boxes.
      //
      {
        for (int i = 0; i < num_v; ++i) {
          sort(boxes[i].begin(), boxes[i].end());
          boxes[i].erase(unique(boxes[i].begin(), boxes[i].end()),
                         boxes[i].end());
        }

        for (V i = 0; i < num_v; ++i)
          for (V j = i + 1; j < num_v; ++j) {
            size_t cnt = 0;
            for (size_t ci = 0, cj = 0;
                 ci < boxes[i].size() && cj < boxes[j].size();) {
              if (boxes[i][ci] == boxes[j][cj]) {
                cnt++;
                ci++;
                cj++;
              } else if (boxes[i][ci] < boxes[j][cj]) {
                ci++;
              } else {
                cj++;
              }
            }
            if (cnt == boxes[j].size()) {
              boxes[j].clear();
            } else if (cnt == boxes[i].size()) {
              boxes[i].clear();
            }
          }
      }

      //
      // Remove unnecessary nodes.
      //
      {
        vector<vector<V>> containd_box_list(num_v);
        for (V i = 0; i < num_v; ++i)
          for (V v : boxes[i]) containd_box_list[v].push_back(i);

        for (V i = 0; i < num_v; ++i)
          sort(containd_box_list[i].begin(), containd_box_list[i].end());

        for (V i = 0; i < num_v; ++i)
          for (V j = i + 1; j < num_v; ++j) {
            if (containd_box_list[i].size() == 0 ||
                containd_box_list[j].size() == 0)
              continue;
            size_t cnt = 0;
            for (size_t ci = 0, cj = 0; ci < containd_box_list[i].size() &&
                                        cj < containd_box_list[j].size();) {
              if (containd_box_list[i][ci] == containd_box_list[j][cj]) {
                cnt++;
                ci++;
                cj++;
              } else if (containd_box_list[i][ci] < containd_box_list[j][cj]) {
                ci++;
              } else {
                cj++;
              }
            }

            if (cnt == containd_box_list[i].size()) {
              containd_box_list[j].clear();
            } else if (cnt == containd_box_list[j].size()) {
              containd_box_list[i].clear();
            }
          }

        boxes.clear();
        boxes.resize(num_v);
        for (V v = 0; v < num_v; ++v)
          for (V b : containd_box_list[v]) boxes[b].push_back(v);

        // Remove pairs of unnecessary twin boxes
        for (V i = 0; i < num_v; ++i)
          for (V j = i + 1; j < num_v; ++j)
            if (containd_box_list[i].size() == 2 &&
                containd_box_list[j].size() == 2) {
              vector<V> &bi1 = boxes[containd_box_list[i][0]];
              vector<V> &bi2 = boxes[containd_box_list[i][1]];
              vector<V> &bj1 = boxes[containd_box_list[j][0]];
              vector<V> &bj2 = boxes[containd_box_list[j][1]];
              if (bi1.size() != 2 || bi2.size() != 2 || bj1.size() != 2 ||
                  bj2.size() != 2)
                continue;

              V k1 = bi1[0] == i ? bi1[1] : bi1[0];
              V k2 = bi2[0] == i ? bi2[1] : bi2[0];
              V l1 = bj1[0] == j ? bj1[1] : bj1[0];
              V l2 = bj2[0] == j ? bj2[1] : bj2[0];

              if (k1 == l1 && k2 == l2) {
                bi2.clear();
                bj1.clear();
              } else if (k1 == l2 && k2 == l1) {
                bi2.clear();
                bj2.clear();
              }
            }
      }

      //
      // Search for boxes that must be contained in the solution.
      //
      {
        vector<vector<V>> containd_box_list(num_v);
        for (V i = 0; i < num_v; ++i)
          for (V v : boxes[i]) containd_box_list[v].push_back(i);
        for (V v = 0; v < num_v; ++v)
          if (containd_box_list[v].size() == 1) {
            V b = containd_box_list[v][0];
            solution.insert(b);
            vector<V> &covered = boxes[b];
            for (V c : covered) containd_box_list[c].clear();
          }

        boxes.clear();
        boxes.resize(num_v);
        for (V v = 0; v < num_v; ++v)
          for (V b : containd_box_list[v]) boxes[b].push_back(v);
      }
    }
  };

  set<V> solution;

  //
  // Create all possible boxes.
  //
  V num_v = g.num_vertices();
  vector<vector<V>> boxes(num_v);
  {
    for (V i = 0; i < num_v; ++i) {
      vector<bool> vis(num_v);
      queue<pair<V, W>> que;

      que.push(make_pair(i, 0));
      vis[i] = true;
      boxes[i].push_back(i);

      while (!que.empty()) {
        V q = que.front().first;
        W dist = que.front().second;
        que.pop();
        if (dist == radius) break;
        for (V v : g.neighbors(q)) {
          if (vis[v]) continue;
          que.push(make_pair(v, dist + 1));
          vis[v] = true;
          boxes[i].push_back(v);
        }
      }
      sort(boxes[i].begin(), boxes[i].end());
    }
  }

  //
  // System split.
  // Find the best solution by DFS
  // The box-covering for tree networks could be performed in O(N^3)
  // while for regular networks it requires O(2^N).
  //
  deque<pair<vector<vector<V>>, set<V>>> que;
  que.push_front(make_pair(boxes, solution));
  size_t min_size = num_v;
  while (!que.empty()) {
    vector<vector<V>> qboxes = que.front().first;
    set<V> qsolution = que.front().second;
    que.pop_front();
    if (qsolution.size() >= min_size) continue;

    burning_splitted(g, radius, qsolution, qboxes);

    coverage_manager cm(g, radius, 1.0);
    for (const V &qv : qsolution) cm.add(g, qv);
    if (cm.get_current_coverage() == 1.0) {
      if (min_size > qsolution.size()) {
        solution = qsolution;
        min_size = solution.size();
      }
      continue;
    }

    // Find the node that is in the smallest number of boxes
    vector<set<V>> containd_box_list(num_v);
    vector<V> covered_largest(num_v, 0);
    for (V i = 0; i < num_v; ++i)
      for (V v : qboxes[i]) {
        containd_box_list[v].insert(i);
        if ((size_t)covered_largest[v] < qboxes[i].size())
          covered_largest[v] = qboxes[i].size();
      }
    V selected_v = -1;
    size_t min_list_size = num_v;
    V max_covered = 0;
    for (int v = 0; v < num_v; ++v) {
      if (containd_box_list[v].empty()) continue;
      if (containd_box_list[v].size() < min_list_size) {
        min_list_size = containd_box_list[v].size();
        max_covered = covered_largest[v];
        selected_v = v;
      } else if (containd_box_list[v].size() == min_list_size &&
                 covered_largest[v] > max_covered) {
        min_list_size = containd_box_list[v].size();
        max_covered = covered_largest[v];
        selected_v = v;
      }
    }

    assert(selected_v >= 0 && containd_box_list[selected_v].size() >= 2);

    // After Selected
    for (V onlycontain : containd_box_list[selected_v]) {
      // Generate subboxes
      vector<vector<V>> subboxes(num_v);
      for (V v = 0; v < num_v; v++) {
        if (v != onlycontain &&
            containd_box_list[selected_v].find(v) !=
                containd_box_list[selected_v].end())
          continue;
        subboxes[v].assign(qboxes[v].begin(), qboxes[v].end());
      }
      que.push_front(make_pair(subboxes, qsolution));
    }
  }
  vector<V> ret(solution.begin(), solution.end());
  return ret;
}

/**
 * Box-Covering algorithm, named Greedy-Coloring.
 * \param g is graph to cover
 * \param diameter is diameter of each box
 */
vector<pair<W, size_t>> box_cover_coloring(const G &g, W diameter) {
  V N = g.num_vertices();
  vector<vector<V>> colors(N, vector<V>(diameter + 1, 0));

  for (V i = 1; i < N; ++i) {  //(iii)
    vector<W> dist(N, N);      //(a)
    queue<pair<V, W>> que;
    que.push({i, 0});
    dist[i] = 0;
    while (!que.empty()) {
      V v = que.front().first;
      W d = que.front().second;
      que.pop();
      for (const auto &u : g.neighbors(v)) {
        if (dist[u] <= d + 1) continue;
        dist[u] = d + 1;
        que.push({u, d + 1});
      }
    }

    multimap<W, V> dict;
    for (V v = 0; v < N; ++v) dict.insert({dist[v], v});

    for (W l_B = 1; l_B <= diameter; ++l_B) {  //(b)
      vector<bool> color_used(N, false);
      for (auto it = dict.lower_bound(l_B); it != dict.end(); it++) {
        V j = it->second;
        if (j >= i) continue;
        color_used[colors[j][l_B]] = true;
      }
      for (int c = 0; c < N; ++c)
        if (!color_used[c]) {
          colors[i][l_B] = c;
          break;
        }
    }
  }

  vector<pair<W, size_t>> ret;
  for (W l_B = 1; l_B <= diameter; ++l_B) {
    set<V> s;
    for (int i = 0; i < N; ++i) s.insert(colors[i][l_B]);
    ret.push_back({l_B, s.size()});
  }
  return ret;
}

/**
 * Box-Covering algorithm, named comapct box burning (CBB).
 * \param g is graph to cover
 * \param diameter is diameter of each box
 */
vector<V> box_cover_cbb(const G &g, W diameter) {
  auto bfs = [](const G &g, W diameter, V center) -> set<V> {
    set<V> ret;
    queue<pair<V, W>> que;
    que.push({center, 0});
    vector<bool> used(g.num_vertices(), false);
    used[center] = true;
    ret.insert(center);
    while (!que.empty()) {
      V v = que.front().first;
      W d = que.front().second;
      que.pop();
      if (d == diameter) break;
      for (V u : g.neighbors(v)) {
        if (used[u]) continue;
        que.push({u, d + 1});
        used[u] = true;
        ret.insert(u);
      }
    }
    return ret;
  };

  vector<V> centers;
  set<V> uncovered;
  for (int i = 0; i < g.num_vertices(); ++i) uncovered.insert(i);

  while (!uncovered.empty()) {
    V center = -1;
    {
      vector<V> tmp(uncovered.begin(), uncovered.end());
      shuffle(tmp.begin(), tmp.end(),agl::random);
      center = tmp[0];
    }
    assert(center >= 0);

    uncovered.erase(center);
    centers.push_back(center);
    set<V> box = bfs(g, diameter, center);
    vector<V> candidate(box.begin(), box.end());
    for (const auto &v : candidate) {
      if (!box.count(v)) continue;
      set<V> next;
      auto b = bfs(g, diameter, v);
      for (const auto &u : b)
        if (box.count(u)) next.insert(u);
      box.swap(next);
    }
    for (const auto &v : box) uncovered.erase(v);
  }
  return centers;
}

/**
 * Box-Covering algorithm, which is using bottom-k min-hash sketch
 * \param g is graph to cover
 * \param radius is radius of each box
 * \param k maximum size of each sketch
 * \param pass the maximum number of iterations for multi-pass improvement
 * \param least_coverage when the coverage exceeds least_coverage, the main loop will terminate
 * \param alpha if the index size is smaller than alpha * k * N, MEMB will be applied
 */
vector<V> box_cover_sketch(const G &g, W radius, const int k, const int pass, double least_coverage, double alpha) {
  coverage_manager cm(g, radius, least_coverage);
  return box_cover_sketch(g, radius, k, pass, cm, alpha);
}

/**
 * Box-Covering algorithm, which is using bottom-k min-hash sketch
 * \param g is graph to cover
 * \param radius is radius of each box
 * \param k maximum size of each sketch
 * \param pass the maximum number of iterations for multi-pass improvement
 * \param cm coverage manager
 * \param alpha if the index size is smaller than alpha * k * N, MEMB will be applied
 */
vector<V> box_cover_sketch(const G &g, W radius, const int k, const int pass, coverage_manager &cm, double alpha) {
  assert(k > 0);

  const V num_v = g.num_vertices();
  const size_t index_size_limt = num_v * k * alpha;

  vector<V> centers;
  vector<V> rank(num_v);
  vector<V> inv(num_v);
  for (V i = 0; i < num_v; ++i) inv[i] = i;

  for (int pass_trial = 0; pass_trial < pass; pass_trial++) {
    shuffle(inv.begin(), inv.end(), agl::random);
    for (int i = 0; i < num_v; ++i) rank[inv[i]] = i;

    bool use_memb = true;
    vector<vector<V>> X =
        build_sketch(g, radius, k, rank, inv, cm, use_memb, index_size_limt);

    if (use_memb) {
      select_lazy_greedily(g, X, rank, inv, centers, cm);
    } else {
      select_greedily(g, X, centers, k, cm);
    }

    if (cm.is_covered()) break;
  }
  return centers;
}
}  // namespace agl
