#include "nsp.hpp"
#include <jian/geom.hpp>
#include <jian/pdb.hpp>
#include <jian/scoring/Scoring.hpp>
#include <jian/utils/file.hpp>

namespace jian {
	namespace {
		void score_s(Scoring & scoring, std::string filename) {
			Chain chain;
			chain_read_model(chain, filename);
			scoring.run(chain);
			std::cout
				<< scoring.m_score_dih << "(dih) "
				<< scoring.m_score_dist << "(dist) "
				<< scoring.m_score << "(total)"
				<< std::endl;
		}

		void score_l(Scoring & scoring, std::string filename) {
			EACH_SPLIT_LINE(filename, " ",
				score_s(scoring, F[0]);
			);
		}

		void train_s(Scoring & scoring, std::string filename) {
			Chain chain;
			chain_read_model(chain, filename);
			scoring.train(chain);
		}

		void train_l(Scoring & scoring, std::string filename) {
			EACH_SPLIT_LINE(filename, " ",
				train_s(scoring, F[0]);
			);
		}

		REGISTER_NSP_COMPONENT(score) {
			Scoring scoring;
			scoring.init();
			if (par.has("train")) {
				if (par.has("s)")) {
					train_s(scoring, par.get("s"));
				}
				else if (par.has("l")) {
					train_l(scoring, par.get("l", "list"));
				}
				scoring.m_dist_anal->print_counts();
			}
			else {
				if (par.has("s")) {
					score_s(scoring, par.get("s"));
				}
				else if (par.has("l")) {
					score_l(scoring, par.get("l", "list"));
				}
			}
		}
	}
} // namespace jian
















