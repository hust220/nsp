#include <iostream>
#include <map>
#include <functional>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>
#include <jian/nuc2d/SSTree.hpp>
#include <jian/utils/exception.hpp>
#include <jian/utils/file.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(view_constraints) {
    auto &&chain = read_model_to_chain(par["s"][0]);
    int i, j;
    std::map<std::string, std::function<double(const Residue &, const Residue &)>> methods {
        {"min", [](const Residue &r1, const Residue &r2){
            double min = 999999;
            double d;
            for (auto && atom1 : r1) {
                for (auto && atom2 : r2) {
                    d = geom::distance(atom1, atom2);
                    if (d < min) {
                        min = d;
                    }
                }
            }
            return min;
        }},
        {"C4*", [](const Residue &r1, const Residue &r2){
            double d;
            auto &atom1 = r1["C4*"];
            auto &atom2 = r2["C4*"];
            d = geom::distance(atom1, atom2);
            return d;
        }}
    };
    auto &method = methods[par["method"][0]];
	BEGIN_READ_FILE(par["c"][0], " ") {
        if (F.size() >= 2) {
            i = JN_INT(F[0]) - 1;
            j = JN_INT(F[1]) - 1;
            std::cout << i + 1 << ' ' << j + 1 << ' ' << method(chain[i], chain[j]) << std::endl;
        }
	} END_READ_FILE;
}

} // namespace jian

