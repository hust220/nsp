#include <jian/dca/Dca.hpp>
#include "nsp.hpp"

namespace jian {
	namespace dca {

		REGISTER_NSP_COMPONENT(dca) {
			int n = 1;
			float step = 1;
			float pw = 0.5;
			std::string mol_type = "RNA";
			std::string method = "mf";

			std::string out_file = par.get("o", "out");
			std::string fa_file = par.get("i", "in");

			par.set(method, "m", "method");
			par.set(n, "n");
			par.set(step, "step");
			par.set(mol_type, "t", "type");
			par.set(pw, "w");

			Dca *dca = FacDca::create(method, mol_type, pw);
			if (method == "mp") dca->set_step(step);
			dca->run(fa_file, out_file, n - 1);
			delete dca;
		}

	}
} // namespace jian

