#include "nsp.hpp"
#include <jian/geom.hpp>
#include <jian/pdb.hpp>
#include <jian/scoring/Scoring.hpp>
#include <jian/utils/file.hpp>

namespace jian {
	namespace {
		void score(Scoring & scoring, std::string filename) {
			Chain chain;
			chain_read_model(chain, filename);

			
			scoring.run(chain);
			std::cout << scoring.m_score_dih << "(dih) "
				<< scoring.m_score_dist << "(dist) "
				<< scoring.m_score << "(total)"
				<< std::endl;
		}
		void train(Scoring & scoring, std::string filename) {

		}
		REGISTER_NSP_COMPONENT(score) {
			Scoring scoring;
			if (par.has("train")) {
				train(scoring, par.get("l", "list");
			} else {
				if (par.has("s")) {
					score(scoring, par.get("s"));
				}
				else if (par.has("l")) {
					EACH_SPLIT_LINE(par.get("l"), " ",
						std::cout << "Score " << F[0] << std::endl;
						score(scoring, F[0]);
					);
				}
			}
			


		}
	}
} // namespace jian
















