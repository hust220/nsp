#include "nsp.hpp"
#include "cg.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(cg_to_aa) {
	Str filename = par.get("s");
	SP<CG> cg = CG::fac_t::make("6p");
	Molecule mol;
	mol_read(mol, filename);
	for (Model &model : mol) {
		for (Chain &chain : model) {
			chain = cg->to_aa(chain);
		}
	}
	JN_OUT << mol;
}

}

