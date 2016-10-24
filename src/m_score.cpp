#include "nsp.hpp"
#include <jian/geom.hpp>
#include <jian/pdb.hpp>
#include <jian/scoring/ScoreAa.hpp>
#include <jian/scoring/ScorePsb.hpp>
#include <jian/utils/file.hpp>

namespace jian {
	namespace {
		void sum_counts(std::string filename, int rows, int cols) {
			Veci v;
			int i, j, d;
			std::ifstream ifile;

			v = Veci::Zero(cols);
			FOPEN(ifile, filename);
			for (i = 0; i < rows; i++) {
				for (j = 0; j < cols; j++) {
					ifile >> d;
					v[j] += d;
				}
			}
			FCLOSE(ifile);

			for (i = 0; i < cols; i++) {
				std::cout << v[i] << ' ';
			}
			std::cout << std::endl;
		}

		void score_s(ScoreBase * scoring, std::string filename) {
			Chain chain;
			chain_read_model(chain, filename);
			scoring->run(chain);
			std::cout <<
				"Score of " << filename << ": " << 
				//scoring->m_score_dih << "(dih) " <<
				//scoring->m_score_dist << "(dist) " <<
				scoring->m_score << "(total)" <<
				std::endl;
		}

		void score_l(ScoreBase * scoring, std::string filename) {
			EACH_SPLIT_LINE(filename, " ",
				score_s(scoring, F[0]);
			);
		}

		void train_s(ScoreBase * scoring, std::string filename) {
			Chain chain;

			std::cout << "Train " << filename << " ..." << std::endl;
			chain_read_model(chain, filename);
			scoring->train(chain);
		}

		void train_l(ScoreBase * scoring, std::string filename) {
			EACH_SPLIT_LINE(filename, " ",
				train_s(scoring, F[0]);
			);
		}

		REGISTER_NSP_COMPONENT(score) {
			if (par.has("sum_counts")) {
				Par::pars_t & pars = par["sum_counts"];
				std::string filename = pars[0];
				int rows = std::stoi(pars[1]);
				int cols = std::stoi(pars[2]);
				sum_counts(filename, rows, cols);
			}
			else {
				std::string method = "aa";

				par.set(method, "method", "m");

				ScoreBase *scoring = FacScorer::create(method);
				scoring->init();
				if (par.has("train")) {
					if (par.has("s)")) {
						train_s(scoring, par.get("s"));
					}
					else if (par.has("l")) {
						train_l(scoring, par.get("l", "list"));
					}
					if (par.has("out")) {
						std::ofstream ofile;
						FOPEN(ofile, par.get("out"));
						scoring->print_counts(ofile);
						FCLOSE(ofile);
					}
					else {
						scoring->print_counts(std::cout);
					}
				}
				else {
					if (par.has("s")) {
						score_s(scoring, par.get("s"));
					}
					else if (par.has("l")) {
						score_l(scoring, par.get("l", "list"));
					}
				}
				delete scoring;
			}
		}
	}
} // namespace jian
















