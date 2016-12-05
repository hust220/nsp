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

BEGIN_JN
	namespace {

		class ClusterComponent {
		public:
			using method_t = std::function<Mat*(const Model &)>;

			Str cg = "C4'";
			Str list_file;
			Str matrix_file;
			Par::pars_t m_trajs;
			int m_bin = 1;
			Str m_method = "kmeans";
			std::shared_ptr<Par> m_par;

			ClusterComponent(const Par &par) {
				par.set(m_method, "method");
				par.set(cg, "cg");
				par.set(list_file, "l", "list");
				par.set(matrix_file, "m", "mat", "matrix");
				par.setv(m_trajs, "trajs");
				par.set(m_bin, "b", "bin");
				m_par = std::make_shared<Par>(par);
			}

			template<typename T, typename U>
			void print_clusters(T &&cluster, U &&names) {
				int n = 0;
				for (auto && clu : cluster.m_clusters) {
					JN_OUT << "Cluster " << n + 1 << " (size: " << clu.size() << "): ";
					for (auto && i : clu) {
						JN_OUT << names[i] << ' ';
					}
					JN_OUT << std::endl;
					n++;
				}
			}

			Mat *read_mat(Str filename) {
					std::ifstream ifile;
					FOPEN(ifile, filename);
					int rows, cols;
					ifile >> rows >> cols;
					Mat *mat = new Mat(rows, cols);
					for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) {
						if (ifile >> (*mat)(i, j)) continue; else throw "jian::mat_from_file failed!";
					}
					FCLOSE(ifile);
					return mat;
			}

			void run() {
				std::map<Str, method_t> methods;
				methods["C4'"] = [](const Model &model) {
							int len = num_residues(model);
							Mat *mat = new Mat(len, 3);
							int n_res = 0;
							for (auto && chain : model) for (auto && res : chain) {
								auto &&atom = res["C4*"];
								for (int i = 0; i < 3; i++) {
									(*mat)(n_res, i) = Num(atom[i]);
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
					Mat *mat = new Mat(len, 3);
					int i = 0;
					for (auto && chain : model) {
						for (auto && res : chain) {
							for (auto && atom : res) {
								if (std::find(atom.name.begin(), atom.name.end(), 'P') == atom.name.end()) {
									for (int k = 0; k < 3; k++) {
										(*mat)(i, k) = Num(atom[k]);
									}
									i++;
								}
							}
						}
					}
					return mat;
				};

				method_t &method = methods[cg];
				std::deque<Mat *> mats;
				Mat *mat;
				std::deque<Str> names;
				auto dist = [](Mat *m1, Mat *m2) {return geom::rmsd(*m1, *m2); };
				auto cluster = Cluster::fac_t::make(m_method, *m_par);
				//Cluster cluster(k);

				if (!list_file.empty()) {
					BEGIN_READ_FILE(list_file, " ") {
						LOG << "Reading: " << F[0] << std::endl;
						auto && model = mol_read_to<Model>(F[0]);
						mats.push_back(method(model));
						names.push_back(model.name);
					} END_READ_FILE;
					LOG << "Clustering..." << std::endl;
					mat = Cluster::to_mat(mats.begin(), mats.end(), dist);
				}
				else if (!m_trajs.empty()) {
                    for (auto && traj : m_trajs) {
                        LOG << "Loading file: " << traj << "..." << std::endl;
                        for_each_model(traj, [this, &names, &mats, &method, &traj](const Model &model, int i) {
                            if (i % m_bin == 0) {
                                LOG << "Reading: model " << i + 1 << " (" << num_residues(model) << " nt)" << std::endl;
                                mats.push_back(method(model));
                                names.push_back(to_str(traj, ":model-", i+1));
                            }
                        });
                    }
					LOG << "Clustering..." << std::endl;
					mat = Cluster::to_mat(mats.begin(), mats.end(), dist);
				}
				else if (!matrix_file.empty()) {
					LOG << "Reading matrix..." << std::endl;
					mat = read_mat(matrix_file);
					for (int i = 0; i < mat->rows(); i++) names.push_back(JN_STR(i+1));
					LOG << "Clustering..." << std::endl;
				}

				(*cluster)(*mat);
				print_clusters(*cluster, names);
				LOG << "All done." << std::endl;

				for (auto && i : mats) delete i;
				delete mat;
			}

		};

		REGISTER_NSP_COMPONENT(cluster) {
			ClusterComponent cluster(par);
			cluster.run();
		}

	}

END_JN

