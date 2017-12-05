#include "nsp.hpp"
#include "test_match_graph.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(match_graph) {
    Graph g1{{0, 1, 2}, {{0, 1, 1}, {0, 2, 1}, {1, 2, 2}}};
    Graph g2{{0, 1, 2, 3}, {{0, 3, 1}, {3, 2, 1}, {2, 1, 2}, {0, 2, 2}}};
    for (auto && match : match_graph(g1, g2)) {
        JN_OUT << "New Match:" << std::endl;
        for (auto && p : match) {
            JN_OUT << p[0] << '-' << p[1] << ' ';
        }
        JN_OUT << std::endl;
    }

}

END_JN

