#include <iostream>
#include <map>
#include <string>
#include <functional>
#include <deque>
#include <Eigen/Dense>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom/suppos.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/Cluster.hpp>
#include <jian/utils/log.hpp>

namespace jian {

	template<typename T, typename U>
	static void print_clusters(T &&cluster, U &&names) {
		int n = 0;
		for (auto && clu : cluster._clusters) {
			OUT << "Cluster " << n + 1 << " (size: " << clu.size() << "): ";
			for (auto && i : clu) {
				OUT << names[i] << ' ';
			}
			OUT << std::endl;
			n++;
		}
	}

	template<typename T>
	static void print_clusters(T &&cluster) {
		int n = 0;
		for (auto && clu : cluster._clusters) {
			OUT << "Cluster " << n + 1 << " (size: " << clu.size() << "): ";
			for (auto && i : clu) {
				OUT << i + 1 << ' ';
			}
			OUT << std::endl;
			n++;
		}
	}

	static Mat * model_to_matrix(const Model &model) {
		int len = num_residues(model);
		Mat *mat = new Mat(len, 3);
		int n_res = 0;
		for (auto && chain : model) for (auto && res : chain) {
			auto &&atom = res["C4*"];
			for (int i = 0; i < 3; i++) {
				(*mat)(n_res, i) = atom[i];
			}
			n_res++;
		}
		return mat;
	}

	static Mat * model_to_matrix_aa(const Model &model) {
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
		Mat *mat = new Mat(len, 3);
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

	static double dist(Eigen::MatrixXd *m1, Eigen::MatrixXd *m2) {
		auto sp = geom::suppos(*m1, *m2);
		return sp.rmsd;
	}

	REGISTER_NSP_COMPONENT(cluster) {
		int k = 5;
		std::string method_name = "cg";
		std::string list_file;
		std::string matrix_file;
		std::string models_file;

		using method_t = std::function<Mat*(const Model &)>;
		static std::map<std::string, method_t> methods{
			{"cg", model_to_matrix},
			{"aa", model_to_matrix_aa}
		};

		par.set(k, "k");
		par.set(method_name, "method");
		par.set(list_file, "l", "list");
		par.set(matrix_file, "m", "mat", "matrix");
		par.set(models_file, "models");

		method_t &method = methods[method_name];
		std::deque<Mat *> mats;
		std::deque<std::string> names;

		Cluster cluster(k);
		if (par.has("list")) {
			BEGIN_READ_FILE(list_file, " ") {
				LOG << "Reading: " << F[0] << std::endl;
				auto && model = mol_read_to<Model>(F[0]);
				mats.push_back(method(model));
				names.push_back(model.name);
			} END_READ_FILE;
			LOG << "Clustering..." << std::endl;
			cluster(mats.begin(), mats.end(), dist);
			print_clusters(cluster, names);
			LOG << "All done." << std::endl;
		}
		else if (par.has("models")) {
			for_each_model(models_file, [&mats, &method](const Model &model, int i) {
				LOG << "Reading: model " << i+1 << " (" << num_residues(model) << " nt)" << std::endl;
				mats.push_back(method(model));
			});
			LOG << "Clustering..." << std::endl;
			cluster(mats.begin(), mats.end(), dist);
			print_clusters(cluster);
			LOG << "All done." << std::endl;
		}
		else if (par.has("matrix")) {
			Mat && mat = mat_from_file(matrix_file);
			cluster(mat);
			print_clusters(cluster);
		}

		for (auto && i : mats) delete i;
	}

} // namespace jian

