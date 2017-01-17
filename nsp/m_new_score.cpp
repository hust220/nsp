#include "nsp.hpp"
#include <jnbio/scoring/new_score.hpp>
#include <jian/utils/file.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(new_score) {
	if (par.has("l")) {
		for (auto &&it : FileLines(par.get("l"))) {
			JN_OUT << scoring::new_score(mol_read_to<Model>(it.arr[0])) << std::endl;
		}
    } else {
		JN_OUT << scoring::new_score(mol_read_to<Model>(par.get("s"))) << std::endl;
    }
}

END_JN
















