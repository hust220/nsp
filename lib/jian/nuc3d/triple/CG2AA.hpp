#pragma once

#include <string>
#include <deque>
#include <vector>
#include <fstream>
#include <map>
#include "CG.hpp"
#include "../../utils/Env.hpp"
#include "../../utils/file.hpp"
#include "../../pdb/Chain.hpp"
#include "../../geom.hpp"

namespace jian {
namespace triple {

class CG2AA {
public:
    using mats_t = std::deque<Mat *>;
    using mats_map_t = std::map<std::string, mats_t>;
    using names_t = std::deque<std::string>;

    mats_map_t d_mats_map;
    names_t d_names;
    std::string d_path;

    CG2AA() {
        d_path = Env::lib() + "/RNA/pars/nuc3d/triple/CG2AA/";
        std::string list_name = d_path + "list";
        EACH_SPLIT_LINE(list_name.c_str(), " ",
            d_mats_map[F[0]].push_back(mat_from_pdb(F[1]));
            d_names.push_back(F[1]);
        );
    }

    ~CG2AA() {
        for (auto && mats : d_mats_map) {
            for (auto && mat : mats.second) {
                delete mat;
            }
        }
    }

    template<typename U>
    Chain operator ()(const Mat & mat, U && seq) {
        Chain chain;
        int n = mat.rows() / 5;
        for (int i = 0; i < n; i++) {
            chain.push_back(make_residue(mat, seq, i));
        }
        return chain;
    }

    Mat *mat_from_pdb(const std::string &name) {
        std::string path = d_path + name;
        Residue r = residue_from_file(path);
        CG::atoms_t &atoms = CG::res_atoms_map[r.name];
        Mat *mat = new Mat(5, 3);
        int i = 0;
        for (auto && atom_name : atoms) {
            for (auto && atom : r) {
                if (atom_name == atom.name) {
                    for (int k = 0; k < 3; k++) {
                        (*mat)(i, k) = atom[k];
                    }
                    i++;
                }
            }
        }
        return mat;
    }

    template<typename U>
    Residue make_residue(const Mat & mat, U && seq, int n) {
        std::string name = seq[n];
        Mat m = mat.block<5, 3>(n * 5, 0);
        mats_t &mats = d_mats_map[name];
        std::vector<double> rmsds(mats.size());
        std::transform(mats.begin(), mats.end(), rmsds.begin(), [&](auto &&m_){
            return geom::rmsd(m, *m_);
        });
        std::vector<double>::iterator it = std::min_element(rmsds.begin(), rmsds.end());
        int d = std::distance(rmsds.begin(), it);
        return residue_from_file(d_path + d_names[d]);
    }

};

template<typename U>
Chain cg_to_aa(const Mat & mat, U && seq) {
    static CG2AA cg2aa;
    return cg2aa(mat, seq);
}

} // namespace triple
} // namespace jian

