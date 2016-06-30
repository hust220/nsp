#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include "CGpsb.hpp"

namespace jian {

int CGpsb::size_res = 3;

Residue CGpsb::res(const Residue &r) {
    static std::vector<std::string> v {"N1", "C2", "N3", "C4", "C5", "C6"};
    Residue res;
    res.name = r.name;
    std::array<double, 3> arr{0, 0, 0};
    for (auto && atom : r) {
        if (atom.name == "C4*") {
            res.push_back(Atom("P", atom[0], atom[1], atom[2]));
        } else if (atom.name == "C1*") {
            res.push_back(Atom("S", atom[0], atom[1], atom[2]));
        } else if (std::find(v.begin(), v.end(), atom.name) != v.end()) {
            for (int i = 0; i < 3; i++) {
                arr[i] += atom[i];
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        arr[i] /= 6.0;
    }
    res.push_back(Atom("B", arr[0], arr[1], arr[2]));
    return res;
}

Chain CGpsb::chain(const Chain &chain) {
    Chain c;
    c.name = chain.name;
    c.model_name = chain.model_name;
    for (auto && r : chain) {
        c.push_back(res(r));
    }
    return c;
}

} // namespace jian
