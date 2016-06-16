#include "nsp.hpp"
#include <jian/pdb/Model.hpp>
#include <jian/geom.hpp>

namespace jian {

static double min_distance(const Residue &r1, const Residue &r2) {
    double d, min = 99999;
    for (auto && a1 : r1) {
        for (auto && a2 : r2) {
            d = geom::distance(a1, a2);
            if (d < min) {
                min = d;
            }
        }
    }
    return min;
}

REGISTER_NSP_COMPONENT(contacts) {
    Chain c = residues_from_file(par[std::vector<std::string>{"pdb", "p"}][0]);
    int len = c.size();
    for (int i = 0; i < len; i++) {
        for (int j = i + 1; j < len; j++) {
            if (min_distance(c[i], c[j]) < std::stoi(par[std::vector<std::string>{"cutoff", "c"}][0])) {
                std::cout << i+1 << ' ' << j+1 << std::endl;
            }
        }
    }
}

} // namespace jian
















