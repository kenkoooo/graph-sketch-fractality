#include "graph.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;
using testing::Types;


namespace {
typedef Types<unweighted_graph> GraphTestTypes;

template<typename T>
class graph_io_test : public testing::Test {};
TYPED_TEST_CASE(graph_io_test, GraphTestTypes);

template<typename GraphType>
void is_graph_eq(const GraphType& g1,const GraphType& g2) {
  EXPECT_TRUE(g1.num_vertices() == g2.num_vertices());
  for(auto v : g1.vertices()) {
    EXPECT_TRUE(g1.degree(v) == g2.degree(v));
    for(size_t i = 0; i < g1.degree(v); i++) {
      const auto& e1 = g1.edge(v,i);
      const auto& e2 = g2.edge(v,i);
      EXPECT_TRUE(to(e1) == to(e2));
      EXPECT_TRUE(is_eq(weight(e1), weight(e2)));
    }
  }
}

TEST(graph_io_test, graph_type_string) {
  EXPECT_TRUE(graph_binary_format<G>() == "unweighted");
}

TYPED_TEST(graph_io_test, binary) {
  using GraphType = TypeParam;
  GraphType g1(add_random_weight<GraphType>(generate_erdos_renyi(100, 2)));
  pretty_print(g1);

  std::ostringstream oss (std::stringstream::binary);
  write_graph_binary(g1, oss);

  std::istringstream iss(oss.str(), std::stringstream::binary);
  auto g2 = read_graph_binary<GraphType>(iss);

  is_graph_eq(g1, g2);
}

TYPED_TEST(graph_io_test, tsv) {
  using GraphType = TypeParam;
  GraphType g1(add_random_weight<GraphType>(generate_erdos_renyi(100, 2)));
  pretty_print(g1);

  std::ostringstream oss (std::stringstream::binary);
  oss << std::setprecision(12);
  write_graph_tsv(g1, oss);

  std::istringstream iss(oss.str(), std::stringstream::binary);
  auto g2 = read_graph_tsv<GraphType>(iss);

  is_graph_eq(g1, g2);
}

TEST(graph_io_test, tsv_failed) {
  string s;
  s += "0 1";
  s += "\n3"; // invalid line(only from vertex)
  s += "\n1 2";
  std::istringstream iss(s);
  ASSERT_DEATH(read_graph_tsv<G>(iss), "CHECK Failed .* at line 2: \"3\".*");
}
} // namespace
