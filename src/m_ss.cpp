#include "nsp.hpp"
#include <jian/scoring/ParBp.hpp>
#include <jian/dca.hpp>

BEGIN_JN
	namespace {
		Str get_ss(const Model &model) {
			auto residues = model.residues();
			ParBp pb;
			int i, j, l;
			::jian::dca::pairs_t pairs;
			
			l = size(residues);
			for (i = 0; i < l; i++) {
				for (j = i + 1; j < l; j++) {
					pb.anal(residues[i], residues[j]);
					if (pb.is_wc()) {
						pairs.push_back({ i, j });
					}
				}
			}
			//::jian::dca::print_pairs(pairs);
			return ::jian::dca::pairs_to_ss(pairs, residues.size());
		}

		REGISTER_NSP_COMPONENT(ss) {
			Par::pars_t mols;

			par.setv(mols, "mols");
			for (auto && mol : mols) {
				for_each_model(mol, [&mol](const Model &model, int i) {
					JN_OUT << to_str(mol, ":model-", i + 1) << std::endl;
					JN_OUT << get_ss(model) << std::endl;
				});
			}
		}
	}
END_JN
















