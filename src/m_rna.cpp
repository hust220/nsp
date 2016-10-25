#include "nsp.hpp"
#include <jian/pdb.hpp>

namespace jian {

	REGISTER_NSP_COMPONENT(rna) {
		std::string in = par.get("s", "pdb", "i", "in");
		std::string out = par.get("o", "out");
		std::ofstream ofile;

		FOPEN(ofile, out);
		ofile << mol_read_to<Model>(in, "RNA") << std::endl;
		FCLOSE(ofile);
	}

	REGISTER_NSP_COMPONENT(dna) {
		std::string in = par.get("s", "pdb", "i", "in");
		std::string out = par.get("o", "out");
		std::ofstream ofile;

		FOPEN(ofile, out);
		ofile << mol_read_to<Model>(in, "DNA") << std::endl;
		FCLOSE(ofile);
	}

} // namespace jian

