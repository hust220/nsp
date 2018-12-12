#include "nsp.hpp"
#include "pdb.hpp"
#include "log.hpp"

namespace jian {
	namespace {

		REGISTER_NSP_COMPONENT(pdb2jn) {
			S in;
			S out;

			par.set(in, "i", "s");
			par.set(out, "o", "out");

			MolSerial(out, std::ios::out).write_mol(in);
		}

		REGISTER_NSP_COMPONENT(jn2pdb) {
			S in = par.get("i", "s");

			MolSerial serial(in, std::ios::in);

			Molecule mol;
			serial.read(mol);
			JN_OUT << mol;
		}

	}
}

