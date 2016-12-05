#pragma once

#include "../Chain.hpp"
#include "../../utils/Cluster.hpp"

BEGIN_JN
namespace pdb {

Cluster::clusters_t cluster_chains(const std::deque<Chain> &chains, int k);

} // namespace pdb
END_JN
