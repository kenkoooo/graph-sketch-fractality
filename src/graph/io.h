#pragma once
#include "base/base.h"
#include "base/io.h"
#include "graph.h"

namespace {
//set the end line code(regardless of the os)
const char kEndLine = '\n';
const std::string kMagic = "AGL_BINARY";
const std::string kVersion = "0.01";

// http://d.hatena.ne.jp/osyo-manga/20120211/1328922379
extern void* enabler;
template<bool B, typename T = void>
using enabler_if = typename std::enable_if<B, T>::type*&;
}

namespace agl {
template<typename GraphType = G>
typename GraphType::edge_list_type read_edge_list_tsv(std::istream &is = std::cin) {
  using E = typename GraphType::E;

  typename GraphType::edge_list_type es;
  std::string line;
  int line_no = 1; //1-origin
  while(std::getline(is,line)) {
    std::istringstream iss(line);
    V v;
    E e;
    CHECK_MSG(iss >> v >> e, ("at line " + to_string(line_no) + ": \"" + line + "\"").c_str());
    es.emplace_back(v,e);
    line_no++;
  }
  return es;
}

template<typename GraphType = G>
typename GraphType::edge_list_type read_edge_list_tsv(const char *filename) {
  if (strcmp(filename, "-") == 0) {
    return read_edge_list_tsv<GraphType>(std::cin);
  } else {
    std::ifstream ifs(filename);
    ifs.sync_with_stdio(false);
    CHECK_PERROR(ifs);
    return read_edge_list_tsv<GraphType>(ifs);
  }
}

template<typename GraphType = G>
GraphType read_graph_tsv(std::istream &is = std::cin) {
  return GraphType(read_edge_list_tsv<GraphType>(is));
}

template<typename GraphType = G>
GraphType read_graph_tsv(const char *filename) {
  return GraphType(read_edge_list_tsv<GraphType>(filename));
}

template<typename EdgeType>
void write_graph_tsv_edge(const EdgeType& e, std::ostream &os){
  os << to(e) << " " << weight(e);
}
template<> void write_graph_tsv_edge(const unweighted_edge& e,std::ostream &os);

template<typename GraphType = G>
void write_graph_tsv(const GraphType &g, std::ostream &os = std::cout) {
  for (auto v : g.vertices()) {
    for (auto e : g.edges(v)) {
      os << v << " ";
      write_graph_tsv_edge(e, os);
      os << std::endl;
    }
  }
}

template<typename GraphType = G>
void write_graph_tsv(const GraphType &g, const char *filename) {
  if (strcmp(filename, "-") == 0) {
    write_graph_tsv(g, std::cout);
  } else {
    std::ofstream ofs(filename);
    CHECK_PERROR(ofs);
    write_graph_tsv(g, ofs);
  }
}

template<typename GraphType = G>
void pretty_print(const GraphType &g, std::ostream &os = std::cerr) {
  static constexpr V kLimitNumVertices = 5;
  static constexpr size_t kLimitNumEdges = 10;

  os << "=========" << std::endl;
  os << "  Vertices: " << g.num_vertices() << std::endl;
  os << "  Edges: " << g.num_edges() << std::endl;
  os << "  Type: " << typename_of(g) << std::endl;
  for (D d : directions()) {
    os << "----------" << std::endl;
    for (V v = 0; v < std::min(kLimitNumVertices, g.num_vertices()); ++v) {
      os << "  " << v << (d == kFwd ? " -> " : " <- ");
      for (size_t i = 0; i < std::min(kLimitNumEdges, g.degree(v, d)); ++i) {
        if (i > 0) os << ", ";
        os << g.edge(v, i, d);
      }
      if (kLimitNumEdges < g.degree(v, d)) os << ", ...";
      os << std::endl;
    }
    if (kLimitNumVertices < g.num_vertices()) {
      os << "  ..." << std::endl;
    }
  }
  os << "=========" << std::endl;
}

//format
template<typename WeightType,
  enabler_if<std::is_integral<WeightType>::value> = enabler >
std::string graph_binary_format_weight() {
  return "weight=int,weight_size=" + to_string(sizeof(WeightType));
}
template<typename WeightType,
  enabler_if<std::is_floating_point<WeightType>::value> = enabler >
std::string graph_binary_format_weight() {
  return "weight=float,weight_size=" + to_string(sizeof(WeightType));
}

template<typename GraphType,
  enabler_if<std::is_same<typename std::decay<GraphType>::type,G>::value> = enabler >
const std::string graph_binary_format() {
  return "unweighted";
}

template<typename GraphType,
  enabler_if<!std::is_same<typename std::decay<GraphType>::type,G>::value> = enabler >
const std::string graph_binary_format() {
  using decayed_gragh_type = typename std::decay<GraphType>::type;
  using weight_type = typename decayed_gragh_type::W;

  return graph_binary_format_weight<weight_type>();
}

// write_binary
template<typename EdgeType,
  enabler_if<std::is_pod<EdgeType>::value> = enabler >
void write_edge_binary(std::ostream& os, const EdgeType& e) {
  write_binary(os, e);
}

template<typename GraphType>
void write_graph_binary(const GraphType &g, std::ostream &os = std::cout) {
  //header
  os << kMagic << kEndLine << kVersion << kEndLine << graph_binary_format<GraphType>() << kEndLine;

  //body
  write_binary(os, g.num_vertices());
  write_binary(os, g.num_edges());
  for(V v = 0; v < g.num_vertices(); v++) {
    write_binary(os, g.degree(v));
    for(std::size_t i = 0; i < g.degree(v); i++) {
      write_binary(os, g.edge(v,i));
    }
  }
  os.flush();
}

template<typename GraphType = G>
void write_graph_binary(const GraphType &g, const char *filename) {
  if (strcmp(filename, "-") == 0) {
    write_graph_binary(g, std::cout);
  } else {
    std::ofstream ofs(filename);
    CHECK_PERROR(ofs);
    write_graph_binary(g, ofs);
  }
}

//reader binary
template<typename EdgeType,
  enabler_if<std::is_pod<EdgeType>::value> = enabler >
EdgeType read_edge_binary(std::istream& is) {
  EdgeType e;
  read_binary(is, &e);
  return e;
}

template<typename GraphType = G>
GraphType read_graph_binary(std::istream &is = std::cin) {
  //header
  std::string magic,version,format;
  std::getline(is, magic, kEndLine);
  std::getline(is, version, kEndLine);
  std::getline(is, format, kEndLine);

  CHECK_MSG(magic == kMagic, "Invalid file magic.");
  CHECK_MSG(version == kVersion, "Invalid file version.");
  CHECK_MSG(format == graph_binary_format<GraphType>(), "Invalid file format.");

  //body
  typename GraphType::V num_vertices;
  std::size_t num_edges;
  read_binary(is, &num_vertices);
  read_binary(is, &num_edges);

  std::vector<std::vector<typename GraphType::E>> edges(num_vertices);

  for(V v = 0; v < num_vertices; v++) {
    std::size_t degree;
    read_binary(is, &degree);
    edges[v].reserve(degree);
    for(std::size_t i = 0; i < degree; i++) {
      edges[v].emplace_back(read_edge_binary<typename GraphType::E>(is));
    }
  }

  GraphType deserialized_graph;
  deserialized_graph.assign(std::move(edges));
  return deserialized_graph;
}

template<typename GraphType = G>
GraphType read_graph_binary(const char *filename) {
  if (strcmp(filename, "-") == 0) {
    return read_graph_binary<GraphType>(std::cin);
  } else {
    std::ifstream ifs(filename);
    ifs.sync_with_stdio(false);
    CHECK_PERROR(ifs);
    return read_graph_binary<GraphType>(ifs);
  }
}
}  // namespace agl
