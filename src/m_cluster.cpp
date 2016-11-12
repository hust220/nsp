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
	namespace {

		template<typename NumType>
		class ClusterComponent {
		public:
			using mat_t = Eigen::Matrix<NumType, -1, -1>;
			using method_t = std::function<mat_t*(const Model &)>;

			int k = 5;
			std::string method_name = "cg";
			std::string list_file;
			std::string matrix_file;
			std::string models_file;
			int m_bin = 1;

			ClusterComponent(const Par &par) {
				par.set(k, "k");
				par.set(method_name, "method");
				par.set(list_file, "l", "list");
				par.set(matrix_file, "m", "mat", "matrix");
				par.set(models_file, "models");
				par.set(m_bin, "b", "bin");
			}

			template<typename T, typename U>
			void print_clusters(T &&cluster, U &&names) {
				int n = 0;
				for (auto && clu : cluster.m_clusters) {
					OUT << "Cluster " << n + 1 << " (size: " << clu.size() << "): ";
					for (auto && i : clu) {
						OUT << names[i] << ' ';
					}
					OUT << std::endl;
					n++;
				}
			}

			mat_t *read_mat(std::string filename) {
					std::ifstream ifile;
					FOPEN(ifile, filename);
					int rows, cols;
					ifile >> rows >> cols;
					mat_t *mat = new mat_t(rows, cols);
					for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) {
						if (ifile >> (*mat)(i, j)) continue; else throw "jian::mat_from_file failed!";
					}
					FCLOSE(ifile);
					return mat;
			}

			void run() {
				std::map<std::string, method_t> methods;
				methods["cg"] = [](const Model &model) {
							int len = num_residues(model);
							mat_t *mat = new mat_t(len, 3);
							int n_res = 0;
							for (auto && chain : model) for (auto && res : chain) {
								auto &&atom = res["C4*"];
								for (int i = 0; i < 3; i++) {
									(*mat)(n_res, i) = NumType(atom[i]);
								}
								n_res++;
							}
							return mat;
				};
				methods["aa"] = [](const Model &model) {
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
					mat_t *mat = new mat_t(len, 3);
					int i = 0;
					for (auto && chain : model) {
						for (auto && res : chain) {
							for (auto && atom : res) {
								if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
									for (int k = 0; k < 3; k++) {
										(*mat)(i, k) = NumType(atom[k]);
									}
									i++;
								}
							}
						}
					}
					return mat;
				};

				method_t &method = methods[method_name];
				std::deque<mat_t *> mats;
				mat_t *mat;
				std::deque<std::string> names;
				auto dist = [](mat_t *m1, mat_t *m2) {return geom::rmsd(*m1, *m2); };
				Cluster cluster(k);

				if (!list_file.empty()) {
					BEGIN_READ_FILE(list_file, " ") {
						LOG << "Reading: " << F[0] << std::endl;
						auto && model = mol_read_to<Model>(F[0]);
						mats.push_back(method(model));
						names.push_back(model.name);
					} END_READ_FILE;
					LOG << "Clustering..." << std::endl;
					mat = Cluster::to_mat<NumType>(mats.begin(), mats.end(), dist);
				}
				else if (!models_file.empty()) {
					for_each_model(models_file, [this, &names, &mats, &method](const Model &model, int i) {
						if (i % m_bin == 0) {
							LOG << "Reading: model " << i + 1 << " (" << num_residues(model) << " nt)" << std::endl;
							mats.push_back(method(model));
							std::ostringstream stream;
							stream << "model-" << i + 1;
							names.push_back(stream.str());
						}
					});
					LOG << "Clustering..." << std::endl;
					mat = Cluster::to_mat<NumType>(mats.begin(), mats.end(), dist);
				}
				else if (!matrix_file.empty()) {
					LOG << "Reading matrix..." << std::endl;
					mat = read_mat(matrix_file);
					for (int i = 0; i < mat->rows(); i++) names.push_back(JN_STR(i+1));
					LOG << "Clustering..." << std::endl;
				}

				cluster(*mat);
				print_clusters(cluster, names);
				LOG << "All done." << std::endl;

				for (auto && i : mats) delete i;
				delete mat;
			}

		};

		REGISTER_NSP_COMPONENT(cluster) {
			std::string precision = "double";

			par.set(precision, "p", "precision");

			if (precision == "float") {
				ClusterComponent<float> cluster(par);
				cluster.run();
			}
			else if (precision == "double") {
				ClusterComponent<double> cluster(par);
				cluster.run();
			}
			else {
				throw "Cluster error! Illegal precision '" + precision + "' !";
			}
		}

	}

} // namespace jian

