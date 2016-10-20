#include "../../geom.hpp"
#include "cluster_chains.hpp"

namespace jian {
namespace pdb {

static Mat *chain_to_matrix_aa(const Chain &chain) {
    int len = 0;
    for (auto && res : chain) {
        for (auto && atom : res) {
            if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
                len++;
            }
        }
    }
    Mat *mat = new Eigen::MatrixXd(len, 3);
    int i = 0;
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
    return mat;
}

static double dist(Mat *m1, Mat *m2) {
    auto sp = geom::suppos(*m1, *m2);
    return sp.rmsd;
}

Cluster::result_t cluster_chains(const std::deque<Chain> &chains, int k) {
    Cluster cluster(k);
    std::deque<Eigen::MatrixXd *> mats; 
    for (auto && chain : chains) {
        mats.push_back(chain_to_matrix_aa(chain));
    }
    auto && result = cluster(mats.begin(), mats.end(), dist);
    for (auto && i : mats) delete i;
    return result;
}


} // namespace pdb
} // namespace jian
