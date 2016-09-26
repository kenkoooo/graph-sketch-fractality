#include "io.h"

namespace agl {
template<> void write_graph_tsv_edge(const unweighted_edge& e,std::ostream &os) {
  os << to(e);
}

}