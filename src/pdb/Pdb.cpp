#include "Pdb.h"

namespace jian {

Pdb::Pdb() {
}

Pdb::Pdb(MolFile &pdb_file) {
    name = pdb_file._name;   
    while (!pdb_file.eof()) {
        models.push_back(Model(pdb_file));
    }
}

//Pdb::Pdb(PdbFile &pdb_file) {
//    name = pdb_file._name;   
//    while (!pdb_file.eof()) {
//        models.push_back(Model(pdb_file));
//    }
//}
//
//Pdb::Pdb(Cif &cif) {
//    name = cif._name;   
//    while (!cif.eof()) {
//        models.push_back(Model(cif));
//    }
//}
//
Pdb::Pdb(string file_name) {
    read(file_name);
}

void Pdb::read(string file_name) {
    if (file_name.size() > 4 && file_name.substr(file_name.size() - 4, 4) == ".pdb") {
        read_pdb(file_name);
    } else if (file_name.size() > 4 && file_name.substr(file_name.size() - 4, 4) == ".cif") {
        Cif cif(file_name);
        (*this) = Pdb(cif);
    } else {
        die("Please give me a file ended with '.pdb' or '.cif'!");
    }
}

void Pdb::read_pdb(string pdbfile) {
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
			models.push_back(Model(lines));
			lines.clear();
		} else if (!line.compare(0, 4, "ATOM")) {
			atom_num++;
			lines.push_back(line);
		}
	}
	if (lines.size() != 0) {
		models.push_back(Model(lines));
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

