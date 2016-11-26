#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/nuc3d/Format.hpp>

namespace jian {
	namespace {
		void m_format(const Par &par, std::string mol_type) {
			std::string in = par.get("s", "pdb", "i", "in");
			Format format;
			Molecule mol;

			mol_read(mol, in, mol_type);
			if (par.has("format")) {
				mol = format(mol);
			}
			else {
				format.sort(mol);
			}
			JN_OUT << mol << std::endl;
		}

		REGISTER_NSP_COMPONENT(format) {
			m_format(par, "");
		}

		REGISTER_NSP_COMPONENT(rna) {
			m_format(par, "RNA");
		}

		REGISTER_NSP_COMPONENT(dna) {
			m_format(par, "DNA");
		}
	}
} // namespace jian

