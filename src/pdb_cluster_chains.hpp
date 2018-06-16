#pragma once

#include "pdb_chain.hpp"
#include "cluster.hpp"

namespace jian {
namespace pdb {

Cluster::clusters_t cluster_chains(const std::deque<Chain> &chains, int k);

} // namespace pdb
}
