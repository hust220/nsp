#include <memory>
#include <jian/pdb.hpp>
#include <jian/cg.hpp>
#include <jian/scoring/ParBp.hpp>
#include "nsp.hpp"

namespace jian {
	namespace {
		REGISTER_NSP_COMPONENT(aa) {
			std::string list_file, prefix;
			Chain chain;
			int l, n, i, j;

			list_file = par.get("l", "list");
			prefix = par.get("prefix");

			n = 0;
			ParBp par_bp;
			//EACH_SPLIT_LINE(list_file, " ",
			//	chain_read_model(chain, F[0]);
			//l = chain.size();
			//for (i = 0; i < l; i++) {
			//	for (j = i + 1; j < l; j++) {
			//		if (par_bp.anal(chain[i], chain[j]).is_paired()) {
			//			Chain c;
			//			c.push_back(chain[i]);
			//			c.push_back(chain[j]);
			//			std::ostringstream stream;
			//			stream << prefix << "." << n + 1 << ".pdb";
			//			mol_write(c, stream.str());
			//			std::cout << n + 1 << ' ' << i + 1 << ' ' << j + 1 << ' ' << par_bp.d << ' ' << par_bp.z << ' ' << par_bp.alpha << std::endl;
			//			std::cout << n + 1 << "wc: " << par_bp.is_wc() << ", nwc: " << par_bp.is_nwc() << std::endl;
			//			n++;
			//		}
			//	}
			//}
			//);
			EACH_SPLIT_LINE(list_file, " ",
				chain_read_model(chain, F[0]);
				l = chain.size();
				for (i = 0; i + 1 < l; i++) {
					j = i + 1;
					par_bp.anal(chain[i], chain[j]);
					Chain c;
					c.push_back(chain[i]);
					c.push_back(chain[j]);
					std::ostringstream stream;
					stream << prefix << "." << n + 1 << ".pdb";
					mol_write(c, stream.str());
					std::cout
						<< n + 1 << ' ' << chain[i].name << i + 1 << ' ' << chain[j].name << j + 1 << ' '
						<< par_bp.d << ' ' << par_bp.z << ' ' << par_bp.alpha << ' ' << par_bp.theta << std::endl;
					//std::cout << n + 1 << "wc: " << par_bp.is_wc() << ", nwc: " << par_bp.is_nwc() << std::endl;
					n++;
				}
			);
		}
	}
} // namespace jian

