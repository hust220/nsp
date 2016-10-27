#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/nuc3d/Format.hpp>

namespace jian {
	namespace {
		void m_format(const Par &par, std::string mol_type) {
			std::string in = par.get("s", "pdb", "i", "in");
			std::string out = par.get("o", "out");
			std::ofstream ofile;
			Format format;
			Molecule mol;

			FOPEN(ofile, out);
			mol_read(mol, in, mol_type);
			format.sort(mol);
			ofile << mol << std::endl;
			FCLOSE(ofile);
		}

		REGISTER_NSP_COMPONENT(rna) {
			m_format(par, "RNA");
		}

		REGISTER_NSP_COMPONENT(dna) {
			m_format(par, "DNA");
		}
	}
} // namespace jian

