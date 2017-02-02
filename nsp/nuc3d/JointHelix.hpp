#pragma once

#include "jian/geom.hpp"
#include "../pdb.hpp"

BEGIN_JN

class JointHelix {
public:
    std::vector<std::string> names {"C5*", "O3*", "C1*"};

    template<typename T, typename U, typename V>
    void set_mat(T &&m, U &&ls, V &&chain) {
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

    template<int N, typename T, typename U>
    Chain _joint(T &&chain1, U &&chain2) {
        int len_chain = chain1.size();
        int len_helix = len_chain / N;
        int num_suppos_atoms = N * 2 * names.size();
        std::deque<int> ls1;
        std::deque<int> ls2;
        for (int i = 0; i < N; i++) {
            ls1.push_back(len_helix * (i + 1) - 2);
            ls1.push_back(len_helix * (i + 1) - 1);
        }
        for (int i = 0; i < N; i++) {
            ls2.push_back(len_helix * i);
            ls2.push_back(len_helix * i + 1);
        }
        Eigen::MatrixXd m1(num_suppos_atoms, 3);
        Eigen::MatrixXd m2(num_suppos_atoms, 3);
        set_mat(m1, ls1, chain1);
        set_mat(m2, ls2, chain2);
        auto sp = geom::suppos(m1, m2);
        INIT_SUPPOS(sp);
        for (auto && res : chain1) {
            for (auto && atom : res) {
                APPLY_SUPPOS(atom, sp);
            }
        }
        Chain chain;
        for (int i = 0; i < N; i++) {
            for (int j = i * len_helix; j < (i + 1) * len_helix; j++) {
                chain.push_back(std::move(chain1[j]));
            }
            for (int j = i * len_helix + 2; j < (i + 1) * len_helix; j++) {
                chain.push_back(std::move(chain2[j]));
            }
        }
        return chain;
    }

    template<int N, typename T, typename U>
    Chain joint(T &&chain1, U &&chain2) {
        int len_chain = chain1.size();
        int len_helix = len_chain / N;
        int num_suppos_atoms = N * 2 * names.size();
        std::deque<int> ls1;
        std::deque<int> ls2;
        for (int i = 0; i < N; i++) {
            if (i % 2 == 0) {
                ls1.push_back(len_helix * (i + 1) - 2);
                ls1.push_back(len_helix * (i + 1) - 1);
            } else {
                ls1.push_back(len_helix * i);
                ls1.push_back(len_helix * i + 1);
            }
        }
        for (int i = 0; i < N; i++) {
            if (i % 2 == 0) {
                ls2.push_back(len_helix * i);
                ls2.push_back(len_helix * i + 1);
            } else {
                ls2.push_back(len_helix * (i + 1) - 2);
                ls2.push_back(len_helix * (i + 1) - 1);
            }
        }
        Eigen::MatrixXd m1(num_suppos_atoms, 3);
        Eigen::MatrixXd m2(num_suppos_atoms, 3);
        set_mat(m1, ls1, chain1);
        set_mat(m2, ls2, chain2);
        auto sp = geom::suppos(m1, m2);
        INIT_SUPPOS(sp);
        for (auto && res : chain1) {
            for (auto && atom : res) {
                APPLY_SUPPOS(atom, sp);
            }
        }
        Chain chain;
        for (int i = 0; i < N; i++) {
            if (i % 2 == 0) {
                for (int j = i * len_helix; j < (i + 1) * len_helix; j++) {
                    chain.push_back(std::move(chain1[j]));
                }
                for (int j = i * len_helix + 2; j < (i + 1) * len_helix; j++) {
                    chain.push_back(std::move(chain2[j]));
                }
            } else {
                for (int j = i * len_helix; j < (i + 1) * len_helix - 2; j++) {
                    chain.push_back(std::move(chain2[j]));
                }
                for (int j = i * len_helix; j < (i + 1) * len_helix; j++) {
                    chain.push_back(std::move(chain1[j]));
                }
            }
        }
        return chain;
    }

};

END_JN






