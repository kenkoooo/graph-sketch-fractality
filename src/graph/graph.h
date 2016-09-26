#pragma once
#include <cstdint>
#include "base/base.h"

namespace agl {
using V = int32_t;
}

#include "direction.h"
#include "edge_type.h"
#include "basic_graph.h"
#include "unweighted_graph.h"

namespace agl {
using G = unweighted_graph;
using E = unweighted_graph::E;
using W = G::W;
static constexpr W kInfW = std::numeric_limits<W>::max();
using unweighted_edge_list = std::vector<std::pair<V, V>>;
}  // namespace agl

#include "io.h"
#include "weight_type.h"
#include "generator.h"
