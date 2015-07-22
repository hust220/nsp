#include "Pdb.h"

namespace jian {

void Pdb::readPDB(string pdbfile) {
	/* set name */
	if (pdbfile.size() > 4 && (pdbfile.substr(pdbfile.size() - 4, 4) == ".pdb" || pdbfile.substr(pdbfile.size() - 4, 4) == ".ent")) {
		name = pdbfile.substr(0, pdbfile.size() - 4);
	} else {
		cerr << "Please give me a file ended with '.pdb' or '.ent'.\n";
		exit(1);
	}

	/* read file */
	ifstream ifile(pdbfile.c_str());
	assert(ifile);

	/* set models */
	int atom_num = 0;
	string line;
	vector<string> lines;
	while (getline(ifile, line, '\n')) {
		if (!line.compare(0, 5, "MODEL")) {
		} else if (!line.compare(0, 6, "ENDMDL")) {
			push(Model(lines));
			lines.clear();
		} else if (!line.compare(0, 4, "ATOM")) {
			atom_num++;
			lines.push_back(line);
		}
	}
	if (lines.size() != 0) {
		push(Model(lines));
		lines.clear();
	}

	if (atom_num == 0) {
		cerr << "PDB '" << name << "' has nothing!" << endl;
		exit(1);
	}

	ifile.close();
}

Model &Pdb::operator [](int n) {
	return models[n];
}

const Model &Pdb::operator [](int n) const {
	return models[n];
}

} /// namespace jian

