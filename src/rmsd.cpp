#include <list>
#include <array>
#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom.hpp>
#include <jian/utils/string.hpp>

namespace jian {

namespace rmsd_detail {

using nums_t = std::vector<int>;

void set_nums(nums_t &nums, const std::string &par) {
    std::vector<std::string> v, w;
    jian::tokenize(par, v, "+");
    for (auto && i : v) {
        jian::tokenize(i, w, "-");
        if (w.size() == 1) {
            nums.push_back(std::stoi(w[0])-1);
        } else if (w.size() == 2) {
            for (int j = std::stoi(w[0]); j <= std::stoi(w[1]); j++) {
                nums.push_back(j-1);
            }
        }
    }
}

} // namespace rmsd_detail

REGISTER_NSP_COMPONENT(rmsd) {
    Chain c1 = residues_from_file(par["pdb"][0]);
    Chain c2 = residues_from_file(par["pdb"][1]);
    assert(c1.size() == c2.size());

    rmsd_detail::nums_t nums;
    if (par.has("nums")) {
        rmsd_detail::set_nums(nums, par["nums"][0]);
    }

    int len = c1.size();
    std::list<std::array<double, 3>> l1, l2;
    int n = 0;
    for (int i = 0; i < len; i++) {
        if (!par.has("nums") || std::find(nums.begin(), nums.end(), i) != nums.end()) {
            if (!par.has("loose")) assert(c1[i].name == c2[i].name);
            for (auto && a1 : c1[i]) {
                for (auto && a2 : c2[i]) {
                    if (a1.name == a2.name) {
                        l1.push_back({a1[0],a1[1],a1[2]});
                        l2.push_back({a2[0],a2[1],a2[2]});
                        n++;
                    }
                }
            }
        }
    }

    Mat m1(n, 3), m2(n, 3);
    int i = 0;
    auto i1 = l1.begin(), i2 = l2.begin();
    for (; i1 != l1.end() && i2 != l2.end(); i1++, i2++, i++) {
        for (int k = 0; k < 3; k++) {
            m1(i, k) = (*i1)[k];
            m2(i, k) = (*i2)[k];
        }
    }
    std::cout << geom::suppos(m1, m2).rmsd << std::endl;
}

} // namespace jian

