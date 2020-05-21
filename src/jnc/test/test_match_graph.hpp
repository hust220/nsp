#include "jian.hpp"

namespace jian {

struct Graph {
    std::vector<int> vertices;
    std::vector<std::vector<int>> edges;
};

using GraphMatch = std::vector<std::array<int, 2>>;

std::vector<GraphMatch> match_graph(const Graph &, const Graph &);

}

