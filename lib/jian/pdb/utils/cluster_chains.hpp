#pragma once

#include "../Chain.hpp"
#include "../../utils/Cluster.hpp"

namespace jian {
namespace pdb {

Cluster::result_t cluster_chains(const std::deque<Chain> &chains, int k);

} // namespace pdb
} // namespace jian
