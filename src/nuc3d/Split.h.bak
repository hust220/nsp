#ifndef SPLIT_H
#define SPLIT_H

#include "../Pdb.h"
#include "loop.h"

namespace jian {

namespace nuc3d {

class Split {
public:
	Split() {
		/* set default parameter */
		char *env_RNA = getenv("RNA");
		lib += env_RNA ? env_RNA : "";
	}
	Split(const Split &split) : lib(split.lib), seq(split.seq), ss(split.ss), 
	                            name(split.name), mol(new RNA(split.mol)), family(split.family) {}
	Split(string parfile) : Split() {
		/* read parameter file */
		ifstream ifile(parfile);
		if (!ifile.is_open()) {
			cerr << "There isn't any file named '" << parfile << "'!" << endl;
			exit(1);
		}
		string str;
		while (ifile >> str) {
			if (!str.compare("pdb_file")) {
				ifile >> str;
				name = str;
			} else if (!str.compare("secondary_structure")) {
				ifile >> str;
				ss = str;
			} else if (!str.compare("library_path")) {
				ifile >> str;
				lib = str + "/";
			} else if (!str.compare("family")) {
				ifile >> str;
				family = str;
			} else {
				cerr << "What does '" << str << "' mean in the parameter file '" << parfile << "'?" << endl;
				exit(1);
			}
		}
		ifile.close();
	
		if (ss == "") {
			cerr << "Please tell me the secondary structure!" << endl;
			exit(1);
		}
		if (name == "") {
			cerr << "Please give me the pdb file!" << endl;
			exit(1);
		}
		mol = new RNA(name);
		seq = mol->getSeq();
	}
	~Split() {
		delete mol;
	}
	Split &operator =(const Split &split) {
		lib = split.lib;
		seq = split.seq;
		ss = split.ss;
		family = split.family;
		name = split.name;
		delete mol;
		mol = new RNA(split.mol);
	}

	void run();
	void splitMol(loop *);
	void writeHelix(helixInfo *);
	void writeLoop(loopInfo *);
	void writeStrand(strandInfo *);

	string lib;
	string ss;
	string family = "other";
	string name;
	
	RNA *mol = NULL;
	string seq;
};

} /// namespace nuc3d

} /// namespace jian


#endif




