#include "Split.h"
#include "Mol2D.h"

void Split::run() {
	Mol2D mol2d(ss, seq);
	splitMol(mol2d.pseudo_head);
}

void Split::splitMol(loop *l) {
	bp *b;
	res *r, *tempRes;
	string outName, tempStr;
	ofstream ofile;
	int i, j, k, t, flag;
	helixInfo *hi;
	loopInfo *li;
	strandInfo *si;
	vector<strand> v;
	int temp1;
	string temp;

	static int iteration = 0;
	iteration++;

	static int n1 = 0;
	static int n2 = 0;
	static int n3 = 0;

	// helix
	if (l->s.head != NULL) {
		n1++;

		stringstream ss;
		ss << "helix_";
		ss << mol->name;
		ss << '_';
		ss << n1;
		ss >> temp;
		outName = lib;
		outName += "helix/";
		outName += temp;
		outName += ".pdb";
		ofile.open(outName.c_str());
		
		hi = new helixInfo;
		hi->name = temp;
		hi->n = l->s.len;
		hi->src = mol->name;

		for (i = 0; i < l->s.getLen(); i++) {
			hi->ss += '(';
		}
		for (i = 0; i < l->s.getLen(); i++) {
			hi->ss += ')';
		}

		for (flag = 0, i = 0; i < (int) mol->chains.size(); i++) {
			for (j = 0; j < (int) mol->chains[i].residues.size(); j++) {
				flag++;
				for (b = l->s.head; b != NULL; b = b->next) {
					if (flag == b->res1.num || flag == b->res2.num) {
						hi->seq += mol->chains[i].residues[j].name;
						for (k = 0; k < (int) mol->chains[i].residues[j].atoms.size(); k++) {
							ofile << mol->chains[i].residues[j].atoms[k].line << endl;
						}
						break;
					}
				}
			}
		}
		ofile.close();
		ss.clear();

		writeHelix(hi);
		//cout << hi->name << ' ' << hi->n << ' ' << hi->seq << ' ' << hi->ss << ' ' << hi->src << endl;
		delete hi;
	}

	// loop
	if (l->head != NULL) {
		n2++;

		stringstream strstr;

		string loop_name;
		strstr << "loop_" << mol->name << '_' << n2;
		strstr >> loop_name;

		//li = new loopInfo(loopName, mol->name, family, l);

		// write to pdb
		strstr.str("");
		strstr.clear();
		string pdb_name;
		strstr << lib << "/loop/" << loop_name << ".pdb";
		strstr >> pdb_name;
		ofile.open(pdb_name.c_str());
		int flag = 0;
		for (int i = 0; i < (int) mol->chains.size(); i++) {
			for (int j = 0; j < (int) mol->chains[i].residues.size(); j++) {
				flag++;
				for (res *r = l->head; r != NULL; r = r->next) {
					if (flag == r->num) {
						for (int k = 0; k < (int) mol->chains[i].residues[j].atoms.size(); k++) {
							ofile << mol->chains[i].residues[j].atoms[k].line << endl;
						}
						break;
					}
				}
			}
		}
		ofile.close();

		// write to info file
		string l_ss = l->getSS();
		int hinge_num = l->hinges.size();
		string info_file;
		strstr.str("");
		strstr.clear();
		if (count_if(l_ss.begin(), l_ss.end(), [](char c) {
			return c == '[' || c == ']';
		})) {
			strstr << lib << "/info/pseudo_knot";
		} else {
			strstr << lib << "/info/loop-" << hinge_num;
		}
		strstr >> info_file;

//		strstr.str("");
//		strstr.clear();
//		string info_file_name;
//		strstr << lib << "/info/loop";
//		strstr >> info_file_name;
		ofile.open(info_file.c_str(), ios::app);
		ofile << loop_name << ' ' << l->getSeq() << ' ' << l->getSS() << ' ' << family << endl;
		ofile.close();

		//writeLoop(li);
		//delete li;
	}
	
	/*
	// strands
	l->getStrand(v);
	for (i = 0; i < (int) v.size(); i++) {
		si = new strandInfo;
		n3++;

		stringstream ss;
		ss << "strand_";
		ss << mol->name;
		ss << '_';
		ss << n3;
		ss >> temp;

		outName = lib;
		outName += "strand/";
		outName += temp;
		outName += ".pdb";
		ss.clear();

		ofile.open(outName.c_str());
		for (flag = 0, j = 0, si->seq = ""; j < (int) mol->chains.size(); j++) {
			for (k = 0; k < (int) mol->chains[j].residues.size(); k++) {
				flag++;
				for (r = v[i].head; r != NULL; r = r->next) {
					if (flag == r->num) {
						si->seq += mol->chains[j].residues[k].name;
						for (t = 0; t < (int) mol->chains[j].residues[k].atoms.size(); t++) {
							ofile << mol->chains[j].residues[k].atoms[t].line << endl;
						}
						break;
					}
				}
			}
		}
		ofile.close();
		
		si->name = temp;
		si->len = v[i].getLen();
		si->dist = v[i].getDist(mol);
		si->src = mol->name;
		writeStrand(si);
		//cout << si->name << ' ' << si->len << ' ' << si->dist  << ' ' << si->seq << ' ' << si->src << endl;
	}
	*/

	// son
	if (l->son != NULL) {
		splitMol(l->son);
	}

	// brother
	if (l->brother != NULL) {
		splitMol(l->brother);
	}
	return;
}

void Split::writeHelix(helixInfo *hi) {
	string name;
	ofstream ofile;

	name = lib + "/info/helix";
	ofile.open(name.c_str(), ios::app);
	ofile << hi->name << ' ' << hi->n << ' ' << hi->seq << ' ' << hi->ss << ' ' << hi->src << endl;
	ofile.close();
}

void Split::writeLoop(loopInfo *li) {
	string name;
	ofstream ofile;

	name = lib + "/info/loop";
	ofile.open(name.c_str(), ios::app);
	ofile << li->name << ' ' << li->n << ' ' << li->flag << ' ' << li->len << ' ' << li->seq << ' ' << li->ss << ' ' << li->src << ' ' << li->family << endl;
	ofile.close();
}

/*
void Split::writeStrand(strandInfo *si) {
	string name;
	ofstream ofile;

	name = lib;
	name += "info/strand";
	ofile.open(name.c_str(), ios::app);
	ofile << si->name << ' ' << si->len << ' ' << si->dist  << ' ' << si->seq << ' ' << si->src << endl;
	ofile.close();
}
*/

