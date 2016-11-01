#include "nsp.hpp"
#include <thread>
#include <jian/nuc3d/Assemble.hpp>
#include <jian/dhmc/DHMC.hpp>
#include <jian/thmc/THMC.hpp>
#include <jian/qhmc/QHMC.hpp>
#include <jian/pdb/utils/cluster_chains.hpp>

namespace jian {
	namespace {
		void jian_3drna_refine(const Par &par, const Chain &chain, int i) {
			try {
				std::string name;
				par.set(name, "name", "job");

				std::ostringstream stream;
				stream << name << ".3drna." << i + 1 << ".log";
				log_file(stream.str());


				nuc3d::DHMC mc;
				mc.init(par);
				mc._pred_chain = chain;

				stream.clear();
				stream.str("");
				stream << name << ".3drna.traj." << i + 1 << ".pdb";
				mc.m_traj = stream.str();

				mc.run();

				stream.clear();
				stream.str("");
				stream << mc._name << ".3drna." << i + 1 << ".pdb";
				mol_write(mc._pred_chain, stream.str());
			}
			catch (const char *c) {
				LOG << c << std::endl;
				exit(1);
			}
		}

		void jian_3drna_tripred(const Par &par, int i) {
			std::string name;
			par.set(name, "name", "job");

			std::ostringstream stream;
			stream << name << ".3drna." << i + 1 << ".log";
			log_file(stream.str());

			nuc3d::triple::THMC tri;
			tri.init(par);
			seed(int(tri._seed + i));

			stream.clear();
			stream.str("");
			stream << name << ".3drna.traj." << i + 1 << ".pdb";
			tri.m_traj = stream.str();

			tri.run();

			stream.clear();
			stream.str("");
			stream << tri._name << ".3drna." << i + 1 << ".pdb";
			mol_write(tri._pred_chain, stream.str());
		}

		REGISTER_NSP_COMPONENT(3drna) {
			int n = 5;
			par.set(n, "n", "num", "number");
			assert(n <= 10);

			std::string ss;
			par.set(ss, "ss", "secondary_structure");

			std::string pred_type;
			par.set(pred_type, "pred_type");

			if (pred_type == "triplex") {
				std::vector<std::thread> t(n);
				//std::thread t[n];
				for (int i = 0; i < n; i++) {
					t[i] = std::thread(jian_3drna_tripred, par, i);
				}
				for (int i = 0; i < n; i++) {
					t[i].join();
				}
			}
			else {
				LOG << "# Initializing Assemble Class..." << std::endl;
				nuc3d::Assemble ass(par);
				LOG << "# Selecting templates..." << std::endl;
				ass.select_templates();

				if (par.has("opt") || std::find_if(ss.begin(), ss.end(), [](auto &&c) {return c != '.' && c != '(' && c != ')'; }) != ss.end() || ass.lack_templates()) {
					std::vector<std::thread> t(n);
					for (int i = 0; i < n; i++) {
						LOG << "# Assembling..." << std::endl;
						ass.assemble();
						LOG << "# Thread " << i + 1 << " for mc..." << std::endl;
						t[i] = std::thread(jian_3drna_refine, par, ass._pred_chain, i);
						LOG << "# Sample all templates..." << std::endl;
						ass.sample_all_templates();
					}
					LOG << "All threads created..." << std::endl;

					for (int i = 0; i < n; i++) {
						t[i].join();
					}
					LOG << "All done..." << std::endl;
				}
				else {
					LOG << "Start sampling..." << std::endl;
					std::deque<Chain> chains;
					for (int i = 0; i < ass._num_sampling; i++) {
						ass.assemble();
						chains.push_back(std::move(ass._pred_chain));
						ass.sample_one_template();
					}
					LOG << "# Totally " << chains.size() << " models sampled." << std::endl;

					LOG << "# Clustering..." << std::endl;
					auto result = pdb::cluster_chains(chains, n);

					for (int i = 0; i < n; i++) {
						std::ostringstream stream;
						stream << ass._name << ".3drna." << i + 1 << ".pdb";
						std::string f = stream.str();
						LOG << "# Writing " << f << "..." << std::endl;
						mol_write(chains[result[i][0]], f);
					}
				}
			}
		}
	}
} // namespace jian

