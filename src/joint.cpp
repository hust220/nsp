#include "nsp.hpp"
#include <jian/nuc3d/JointHelix.hpp>
#include <jian/pdb/Model.hpp>
#include <jian/geom.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(joint) {
    auto set_nums = [](auto &&ls, auto &&v){
        for (auto && s : v) {
            ls.push_back(JN_INT(s)-1);
        }
    };
    auto set_mat = [](auto &&m, auto &&ls, auto &&chain){
        static std::vector<std::string> names {"C5*", "O3*", "C1*"};
        for (int i = 0; i < ls.size(); i++) {
            for (auto && atom : chain[ls[i]]) {
                auto it = std::find(names.begin(), names.end(), atom.name);
                if (it != names.end()) {
                    int j = std::distance(names.begin(), it);
                    for (int k = 0; k < 3; k++) {
                        m(i * 3 + j, k) = atom[k];
                    }
                }
            }
        }
    };
    auto &&chain1 = residues_from_file(par["pdb1"][0]);
    auto &&chain2 = residues_from_file(par["pdb2"][0]);
    std::deque<int> ls1;
    std::deque<int> ls2;
    set_nums(ls1, par["num1"]);
    set_nums(ls2, par["num2"]);
    int len = ls1.size();
    Eigen::MatrixXd m1(len * 3, 3);
    Eigen::MatrixXd m2(len * 3, 3);
    set_mat(m1, ls1, chain1);
    set_mat(m2, ls2, chain2);
    auto sp = geom::suppos(m1, m2);
    INIT_SUPPOS(sp);
    for (auto && res : chain1) {
        for (auto && atom : res) {
            APPLY_SUPPOS(atom, sp);
        }
    }
    for (auto && res : chain2) {
        chain1.push_back(std::move(res));
    }
    residues_to_file(chain1, par["out"][0]);
}

REGISTER_NSP_COMPONENT(joint3) {
    auto &&chain1 = residues_from_file(par[2]);
    auto &&chain2 = residues_from_file(par[3]);
    JointHelix joint_helix;
    std::cout << joint_helix.joint<3>(chain1, chain2) << std::endl;
}

REGISTER_NSP_COMPONENT(joint4) {
    auto &&chain1 = residues_from_file(par[2]);
    auto &&chain2 = residues_from_file(par[3]);
    JointHelix joint_helix;
    std::cout << joint_helix._joint<4>(chain1, chain2) << std::endl;
}

} // namespace jian
















