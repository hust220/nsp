#ifndef R5P_H
#define R5P_H

#include "RNA.h"

namespace jian {

class R5P: public RNA {
public:
	R5P() {}
	R5P(const Model &model): RNA(model) {
		init();
	}
//	R5P(const RNA &rna): RNA(rna) {
//		init();
//	}

private:
	void init() {
		for (auto &chain: chains) {
			for (auto &residue: chain.residues) {
				Residue new_res(residue);
				residue.atoms.clear();
				if (residue.name == "A") {
//					Point p1, p2;
					for (auto &atom: new_res.atoms) {
//						if (set<string>{"C5*", "O3*", "C1*", "N6"}.count(atom.name)) {
						if (set<string>{"C5*", "O3*", "C1*", "N6", "C2"}.count(atom.name)) {
							residue.atoms.push_back(atom);
//						} else if (atom.name == "C2") {
//							p1 = Point(atom.x, atom.y, atom.z);
//						} else if (atom.name == "C5") {
//							p2 = Point(atom.x, atom.y, atom.z);
						}
					}
//					Point p3 = Point(
//						1.5 * p1.x - 0.5 * p2.x, 
//						1.5 * p1.y - 0.5 * p2.y, 
//						1.5 * p1.z - 0.5 * p2.z
//					);
//					residue.atoms.push_back(Atom(p3, "O2"));
				} else if (residue.name == "U") {
					for (auto &atom: new_res.atoms) {
						if (set<string>{"C5*", "O3*", "C1*", "O2", "O4"}.count(atom.name)) {
							residue.atoms.push_back(atom);
						}
					}
				} else if (residue.name == "G") {
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

} /// namespace R5P

#endif // R5P_H

