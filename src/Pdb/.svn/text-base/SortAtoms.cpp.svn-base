#include "SortAtoms.h"

namespace jian {

void SortAtoms::operator ()(Residue &res) {
	std::sort(res.atoms.begin(), res.atoms.end(), [&](
	const Atom &atom1, const Atom &atom2) {
		return _atom_rank[res.name][atom1.name] < _atom_rank[res.name][atom2.name];
	});
}

void SortAtoms::operator ()(Chain &chain) {
	for (auto &&residue: chain.residues) {
		(*this)(residue);
	}
}

void SortAtoms::operator ()(Model &mol) {
	for (auto &&chain: mol.chains) {
		(*this)(chain);
	}
}

void SortAtoms::operator ()(Pdb &pdb) {
	for (auto &&mol: pdb.models) {
		(*this)(mol);
	}
}

} /// namespace jian


