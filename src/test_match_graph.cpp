#include "test_match_graph.hpp"

BEGIN_JN

struct GraphMatching {
    std::vector<GraphMatch> matches;

    GraphMatching(const Graph &g1, const Graph &g2) {
        set_conn(g1, g2);
        for (int i = 0; i < g2.vertices.size(); i++) {
            matches.push_back(get_match(g1, g2, i));
        }
    }

private:
    using Conn = std::vector<std::vector<int>>;
    Conn m_conn1;
    Conn m_conn2;

    void malloc_conn(Conn &conn, int n) {
        conn.resize(n);
        for (int i = 0; i < n; i++) conn[i].resize(n, 0);
    }

    void set_conn(Conn &conn, const Graph &g) {
        for (auto && e : g.edges) {
            conn[e[0]][e[1]] = e[2];
            conn[e[1]][e[0]] = e[2];
        }
    }

    void set_conn(const Graph &g1, const Graph &g2) {
        int n1 = g1.vertices.size();
        int n2 = g2.vertices.size();

        malloc_conn(m_conn1, n1);
        malloc_conn(m_conn2, n2);

        set_conn(m_conn1, g1);
        set_conn(m_conn2, g2);
    }

    GraphMatch get_match(const Graph &g1, const Graph &g2, int i) {
        GraphMatch match;
        bool flag;

        int n1 = g1.vertices.size();
        int n2 = g2.vertices.size();

        match.push_back({0, i});

        do {
            flag = false;
            for (int i1 = 0; i1 < n1; i1++) {
                for (int i2 = 0; i2 < n2; i2++) {
                    if (!in_match(match, i1, i2)) {
                        if (match.empty() || is_match(g1, g2, match, i1, i2)) {
                            flag = true;
                            match.push_back({i1, i2});
                        }
                    }
                }
            }
        } while (flag);

        return match;
    }

    bool in_match(const GraphMatch &m, int i1, int i2) {
        return std::find_if(m.begin(), m.end(), [&i1](auto && p){return p[0] == i1;}) != m.end() ||
               std::find_if(m.begin(), m.end(), [&i2](auto && p){return p[1] == i2;}) != m.end();
    }

    bool is_match(const Graph &g1, const Graph &g2, const GraphMatch &m, int i1, int i2) {
        bool flag = false;
        for (auto && p : m) {
            int a = m_conn1[p[0]][i1];
            int b = m_conn2[p[1]][i2];
            if (a != 0 && b != 0) {
                if (a != b) return false;
                else flag = true;
            }
        }
        return flag;
    }
};

std::vector<GraphMatch> match_graph(const Graph &g1, const Graph &g2) {
    return GraphMatching(g1, g2).matches;
}

END_JN

