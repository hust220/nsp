#include "geom.hpp"
#include "pdb_cluster_chains.hpp"

BEGIN_JN

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
        return geom::rmsd(*m1, *m2);
	}

	Cluster::clusters_t cluster_chains(const std::deque<Chain> &chains, int k) {
		auto cluster = Cluster::fac_t::make("kmeans", Par("k", k));
		std::deque<Mat *> mats;
		for (auto && chain : chains) {
			mats.push_back(chain_to_matrix_aa(chain));
		}
		(*cluster)(mats.begin(), mats.end(), dist);
		for (auto && i : mats) delete i;
		return cluster->m_clusters;
	}


} // namespace pdb
END_JN
