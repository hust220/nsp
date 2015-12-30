#ifndef D5P_H
#define D5P_H

#include "DNA.h"

namespace jian {

class D5P: public DNA {
public:
	D5P() {}
	D5P(string str): D5P(Model(str)) {}
	D5P(const Model &model): DNA(model) {
		init();
	}
//	D5P(const DNA &rna): DNA(rna) {
//		init();
//	}

private:
	void init() {
		for (auto &chain: chains) {
			for (auto &residue: chain.residues) {
				Residue new_res(residue);
				residue.atoms.clear();
				if (residue.name == "DA") {
					for (auto &atom: new_res.atoms) {
						if (set<string>{"C5*", "O3*", "C1*", "N6", "C2"}.count(atom.name)) {
							residue.atoms.push_back(atom);
						}
					}
				} else if (residue.name == "DT") {
					for (auto &atom: new_res.atoms) {
						if (set<string>{"C5*", "O3*", "C1*", "O2", "O4"}.count(atom.name)) {
							residue.atoms.push_back(atom);
						}
					}
				} else if (residue.name == "DG") {
					for (auto &atom: new_res.atoms) {
						if (set<string>{"C5*", "O3*", "C1*", "O6", "N2"}.count(atom.name)) {
							residue.atoms.push_back(atom);
						}
					}
				} else {
					for (auto &atom: new_res.atoms) {
						if (set<string>{"C5*", "O3*", "C1*", "O2", "N4"}.count(atom.name)) {
							residue.atoms.push_back(atom);
						}
					}
				}
			}
		}
	}
};

} /// namespace jian

#endif // D5P_H

