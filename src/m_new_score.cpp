#include "nsp.hpp"
#include <jian/scoring/new_score.hpp>
#include <jian/utils/file.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(new_score) {
	if (par.has("l")) {
		BEGIN_READ_FILE(par["l"][0], " ") {
			std::cout << scoring::new_score(mol_read_to<Model>(F[0])) << std::endl;
		} END_READ_FILE;
    } else {
        std::cout << scoring::new_score(mol_read_to<Model>(par["s"][0])) << std::endl;
    }
}

END_JN
















