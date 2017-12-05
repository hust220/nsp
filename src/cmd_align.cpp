#include <iostream>
#include <map>
#include <string>
#include <functional>
#include <deque>
#include <Eigen/Dense>
#include "nsp.hpp"
#include "pdb.hpp"
#include "cluster.hpp"
#include "geom_suppos.hpp"
#include "file.hpp"

BEGIN_JN

namespace {

    Eigen::MatrixXd * model_to_mat_aa(const Model &model) {
        int len = 0;
        for (auto && chain : model) {
            for (auto && res : chain) {
                for (auto && atom : res) {
                    if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
                        len++;
                    }
                }
            }
        }
        Eigen::MatrixXd *mat = new Eigen::MatrixXd(len, 3);
        int i = 0;
        for (auto && chain : model) {
            for (auto && res : chain) {
                for (auto && atom : res) {
                    if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
                        for (int k = 0; k < 3; k++) {
                            (*mat)(i, k) = atom[k];
                        }
                        i++;
                    }
                }
            }
        }
        return mat;
    }

    template<typename T> Vector<Int> get_nums(T && strs) {
        Vector<Int> v;
        Vector<Str> w;
        for (auto && str : strs) {
            tokenize(str, w, "-");
            if (w.size() == 2) {
                for (int i = JN_INT(w[0]); i <= JN_INT(w[1]); i++) {
                    v.push_back(i-1);
                }
            }
            else if (w.size() == 1) {
                v.push_back(JN_INT(w[0])-1);
            }
        }
        return std::move(v);
    }

    template<typename T> Mat get_mat(const Model &m, T && nums) {
        int n = size(nums);
        Mat mat(n * 3, 3);
        auto && rs = m.residues();
        for (int i = 0; i < n; i++) {
            auto & atom1 = rs[nums[i]]["C5*"];
            auto & atom2 = rs[nums[i]]["O3*"];
            auto & atom3 = rs[nums[i]]["C1*"];
            for (int j = 0; j < 3; j++) {
                mat(i * 3, j) = atom1[j];
                mat(i * 3 + 1, j) = atom2[j];
                mat(i * 3 + 2, j) = atom3[j];
            }
        }
        return std::move(mat);
    }

    REGISTER_NSP_COMPONENT(align_ts) {
        auto && model1 = mol_read_to<Model>(par.get("s1"));
        auto && model2 = mol_read_to<Model>(par.get("s2"));
        auto && n1 = get_nums(par.getv("n1"));
        auto && n2 = get_nums(par.getv("n2"));
        auto && mat1 = get_mat(model1, n1);
        auto && mat2 = get_mat(model2, n2);
        auto && sp = geom::suppos(mat1, mat2);
        for (auto && chain1 : model1) {
            for (auto && res1 : chain1) {
                for (auto && atom1 : res1) {
                    sp.apply(atom1);
                }
            }
        }
        JN_OUT << model1 << std::endl;
    }

    REGISTER_NSP_COMPONENT(align) {
        Model m1, m2;
        mol_read(m1, par["s"][0]);
        mol_read(m2, par["s"][1]);
        Eigen::MatrixXd *mat1 = model_to_mat_aa(m1);
        Eigen::MatrixXd *mat2 = model_to_mat_aa(m2);
        auto sp = geom::suppos(*mat1, *mat2);
        INIT_SUPPOS(sp);
        for (auto && chain : m1) {
            for (auto && res : chain) {
                for (auto && atom : res) {
                    APPLY_SUPPOS(atom, sp);
                }
            }
        }
        std::cout << m1 << std::endl;
        delete mat1;
        delete mat2;
    }

}

END_JN

