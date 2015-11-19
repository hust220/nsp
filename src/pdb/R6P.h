#ifndef R6P_H
#define R6P_H

#include "RNA.h"

namespace jian {

class R6P: public RNA {
public:
	R6P() {}
	R6P(const RNA &rna): RNA(rna) {
		init();
	}

private:
	void init() {
		for (auto &chain: chains) {
			for (auto &residue: chain.residues) {
				Residue new_res(residue);
				residue.atoms.clear();
				for (auto &atom: new_res.atoms) {
					if (set<string>{"P", "O5*", "C5*", "C4*", "C3*", "O3*"}.count(atom.name)) {
						residue.atoms.push_back(atom);
					}
				}
			}
		}
	}
};

} /// namespace jian

#endif // R6P_H

