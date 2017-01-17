#include "nsp.hpp"
#include <jnbio/pdb.hpp>
#include <jnbio/nuc3d/Format.hpp>

BEGIN_JN
	namespace {
		void m_format(const Par &par, S mol_type) {
			S in = par.get("s", "pdb", "i", "in");
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
END_JN

