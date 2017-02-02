#include <vector>
#include <map>
#include <string>
#include <memory>
#include "jian/geom.hpp"
#include "jian/utils/Env.hpp"
#include "jian/utils/file.hpp"
#include "new_score.hpp"

BEGIN_JN

class NewScore {
public:
	using Cluster = struct { double f, beta; std::unique_ptr<Eigen::MatrixXd> center; };
	using Clusters = std::vector<Cluster>;
	using AllClusters = std::array<Clusters, 16>;

	static NewScore &instance() {
		static NewScore ns;
		return ns;
	}

	AllClusters d_all_clusters;
	std::map<std::string, int> d_res_map{
		{"A", 0}, {"U", 1}, {"G", 2}, {"C", 3},
		{"DA", 0}, {"DT", 1}, {"DG", 2}, {"DC", 3}
	};
	std::vector<std::vector<std::string>> d_coarse_grained_atoms{
		{"C5*", "O3*", "C1*", "C2", "C6"},
		{"C5*", "O3*", "C1*", "C2", "C4"},
		{"C5*", "O3*", "C1*", "C2", "C6"},
		{"C5*", "O3*", "C1*", "C2", "C4"}
	};


	NewScore() {
		std::vector<std::string> v{ "A", "U", "G", "C" };
		S path = Env::lib() + "/RNA/pars/scoring/new_score/";
		int i = 0;
		int n;
		S str, list_name, pdb_name;
		double sum = 0;
		for (auto && clusters : d_all_clusters) {
			str = v[i / 4] + v[i % 4];
			list_name = path + str + "/list";
			std::deque<int> lens;
			std::deque<std::string> names;
			for (auto &&it : FileLines(list_name)) {
				if (size(it.arr) == 2) {
					n = std::stoi(it.arr[0]);
					sum += n;
					lens.push_back(n);
					names.push_back(it.arr[1]);
				}
			}
			int len = lens.size();
			clusters.resize(len);
			for (int j = 0; j < len; j++) {
				clusters[j].f = lens[j];
				clusters[j].beta = std::log(lens[j]);
				clusters[j].center.reset(mat_from_pdb(path + str + "/" + std::to_string(j + 1) + "/" + names[j]));
			}
			i++;
		}
		for (auto && clusters : d_all_clusters) {
			for (auto && cluster : clusters) {
				cluster.f = cluster.f / sum;
			}
		}
	}

	void set_res_coords(Eigen::MatrixXd &mat, int beg, const Residue &r) {
		auto &v = d_coarse_grained_atoms[d_res_map[r.name]];
		for (auto && atom : r) {
			auto t = std::find(v.begin(), v.end(), atom.name);
			if (t != v.end()) {
				auto d = std::distance(v.begin(), t);
				for (int i = 0; i < 3; i++) {
					mat(beg + d, i) = atom[i];
				}
			}
		}
	}

	Eigen::MatrixXd *mat_from_pdb(const S &s) {
		auto && m = mol_read_to<Model>(s);
		int len = 5 * num_residues(m);
		Eigen::MatrixXd *mat = new Eigen::MatrixXd(len, 3);
		int n_res = 0;
		for (auto && chain : m) {
			for (auto && res : chain) {
				set_res_coords(*mat, n_res * 5, res);
				n_res++;
			}
		}
		return mat;
	}

	int index_bp(const Residue &r1, const Residue &r2) {
		return d_res_map[r1.name] * 4 + d_res_map[r2.name];
	}

	Eigen::MatrixXd *mat_bp(const Residue &r1, const Residue &r2) {
		Eigen::MatrixXd *mat = new Eigen::MatrixXd(10, 3);
		set_res_coords(*mat, 0, r1);
		set_res_coords(*mat, 5, r2);
		return mat;
	}

	double en_bp(const Residue &r1, const Residue &r2) {
		std::unique_ptr<Eigen::MatrixXd> mat{ mat_bp(r1, r2) };
		int index = index_bp(r1, r2);
		return en_bp(*mat, index);
	}

	double en_bp(const Eigen::MatrixXd &mat, int index) {
		double rmsd;
		double en = 0;
		for (auto && cluster : d_all_clusters[index]) {
			rmsd = geom::rmsd(mat, *(cluster.center));
			en += -cluster.f * std::exp(-cluster.beta * rmsd);
		}
		return en;
	}

	Eigen::VectorXd center(const Residue &r) {
		Eigen::VectorXd v = Eigen::VectorXd::Zero(3);
		double n = 0;
		for (auto && atom : r) {
			for (int i = 0; i < 3; i++) {
				v[i] += atom[i];
			}
			n++;
		}
		return v / n;
	}

	double min_distance(const Residue &r1, const Residue &r2) {
		double d, min = 999;
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

	double score(const Model &model) {
		int num_res1, num_res2;
		double en = 0;

		num_res1 = 0;
		for (auto && chain1 : model) {
			for (auto && res1 : chain1) {
				num_res2 = 0;
				for (auto && chain2 : model) {
					for (auto && res2 : chain2) {
						if (num_res2 > num_res1) {
							//                            if (min_distance(res1, res2) < 4) {
							if (geom::distance(center(res1), center(res2)) < 20) {
								en += en_bp(res1, res2);
							}
						}
						num_res2++;
					}
				}
				num_res1++;
			}
		}
		return en;
	}
};

double scoring::new_score(const Model &model) {
	return NewScore::instance().score(model);
}

double scoring::new_score(const Eigen::MatrixXd &mat, int index) {
	return NewScore::instance().en_bp(mat, index);
}

END_JN

