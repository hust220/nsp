#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/utils/log.hpp>

namespace jian {
	namespace {

		REGISTER_NSP_COMPONENT(pdb2jn) {
			std::string in;
			std::string out;

			par.set(in, "i", "s");
			par.set(out, "o", "out");

			MolSerial(out, std::ios::out).write_mol(in);
		}

		REGISTER_NSP_COMPONENT(jn2pdb) {
			std::string in = par.get("i", "s");

			MolSerial serial(in, std::ios::in);

			Molecule mol;
			serial.read(mol);
			JN_OUT << mol;
		}

	}
} // namespace jian

