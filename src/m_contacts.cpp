#include <sstream>
#include "nsp.hpp"
#include <nsp/pdb.hpp>
#include <jian/geom.hpp>

BEGIN_JN

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
    auto g = par.getv("global");
    Str filename = g[1];
    Num cutoff = 3.5;
    par.set(cutoff, "c", "cutoff");

    Chain c = read_model_to_chain(filename);
    int len = c.size();
    for (int i = 0; i < len; i++) {
        for (int j = i + 3; j < len; j++) {
            if (min_distance(c[i], c[j]) < cutoff) {
                JN_OUT << i+1 << ' ' << c[i].name << ' ' << j+1 << ' ' << c[j].name << std::endl;
                if (par.has("write")) {
                    std::ostringstream stream;
                    stream << c.model_name << '-' << i+1 << '-' << j+1 << ".pdb";
                    mol_write(make_chain({c[i], c[j]}), stream.str());
                }
            }
        }
    }
}

END_JN
















