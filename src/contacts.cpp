#include <sstream>
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

template<typename T>
static Chain make_chain(std::initializer_list<T> ls) {
    Chain chain;
    for (auto && r : ls) {
        chain.push_back(r);
    }
    return chain;
}

REGISTER_NSP_COMPONENT(contacts) {
    Chain c = residues_from_file(par[std::vector<std::string>{"pdb", "p"}][0]);
    int len = c.size();
    for (int i = 0; i < len; i++) {
        for (int j = i + 3; j < len; j++) {
            if (min_distance(c[i], c[j]) < std::stoi(par[std::vector<std::string>{"cutoff", "c"}][0])) {
                std::cout << i+1 << ' ' << c[i].name << ' ' << j+1 << ' ' << c[j].name << std::endl;
                if (par.has("write")) {
                    std::ostringstream stream;
                    stream << c.model_name << '-' << i+1 << '-' << j+1 << ".pdb";
                    residues_to_file(make_chain({c[i], c[j]}), stream.str());
                }
            }
        }
    }
}

} // namespace jian
















