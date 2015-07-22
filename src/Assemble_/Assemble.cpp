#include "Assemble.h"

Assemble::Assemble(char *parfile) {
	/* set default parameter */
	num = 10;
	lib = getenv("RNA");
	job = "assemble";
	out = "./";
	logfile = "3dRNA.log";
	results = NULL;
	family = "other";

	/* read parameter file */
	ifstream ifile(parfile);
	if (!ifile.is_open()) {
		cerr << "I have no found any file named '" << parfile << "'!" << endl;
		exit(1);
	}
	string str;
	while (ifile >> str) {
		if (!str.compare("sequence")) {
			ifile >> str;
			seq = str;
			for (int i = 0; i < seq.size(); i++) {
				if (seq[i] > 90) {
					seq[i] -= 32;
				}
			}
		} else if (!str.compare("secondary_structure")) {
			ifile >> str;
			ss = str;
		} else if (!str.compare("library_path")) {
			ifile >> str;
			lib = str;
		} else if (!str.compare("job_name")) {
			ifile >> str;
			job = str;
		} else if (!str.compare("output_directory")) {
			ifile >> str;
			out = str;
		} else if (!str.compare("number")) {
			ifile >> str;
			num = atoi(str.c_str());
		} else if (!str.compare("view")) {
			ifile >> str;
			_view = stoi(str);
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
	if (seq == "") {
		cerr << "Please tell me the sequence!" << endl;
	}

	/* check the length of the secondary structure and sequence */
	int n = 0;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '.' || ss[i] == '(' || ss[i] == ')' || ss[i] == '[' || ss[i] == ']') {
			n++;
		}
	}
	if (n != seq.size()) {
		cerr << "The length of the secondary structure and sequence should be equal!" << endl;
		exit(1);
	}
}

Assemble::Assemble(string ss, string seq) {
	/* set default parameter */
	num = 10;
	lib = getenv("RNA");
	job = "assemble";
	out = "./";
	logfile = "3dRNA.log";
	results = NULL;
	family = "other";

	/* set secondary structure and sequence */
	this->ss = ss;
	this->seq = seq;

	/* check the length of secondary structure and sequence */
	int n = 0;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '.' || ss[i] == '(' || ss[i] == ')' || ss[i] == '[' || ss[i] == ']') {
			n++;
		}
	}
	if (n != seq.size()) {
		cerr << "The length of the secondary structure and sequence should be equal!" << endl;
		exit(1);
	}

}

Assemble::~Assemble() {
	if (mol != NULL) delete mol;
}

string Assemble::output_path() {
	return out;
}

string Assemble::job_name() {
	return job;
}

RNAs *Assemble::getResults() {
	return results;
}

void Assemble::log(string contents) {
	ofstream ofile(logfile.c_str(), ios::app);
	ofile << contents << endl;
	ofile.close();
}

void Assemble::run() {
	/* assemble */
	ofstream ofile(logfile.c_str(), ios::app);
	ofile << "job: " << job << endl;
	ofile << Time::getTime() << "  start..." << endl;
	ofile << seq << endl;
	ofile << ss << endl;
	ofile.close();

	if (mol != NULL) delete mol;
	mol = new Mol2D(ss, seq);

	results = findRNA(mol->pseudo_head, num);
	if (results->getLen() > num) {
		results->resize(num);
	}

	// compute counts of chain
	int temp = 0;
	for (int i = 0; i < (int) ss.size(); i++) {
		if (ss[i] == '&' || (ss[i] == ')' && i - 1 >= 0 && ss[i - 1] == '(')) {
			temp++;
		}
	}
	temp++;

	// construct new chains
	log("  construct new chains...");
	for (int i = 0; i < results->getLen(); i++) {
		RNA *rna = results->RNAList[i];
		Chain *chain = new Chain[temp];
		int i1 = 0;
		for (int j = 0, k = 0; i1 < (int) rna->chains.size(); i1++) {
			for (int i2 = 0; i2 < (int) rna->chains[i1].residues.size();) {
				if (ss[k] == '&' || (ss[k] == ')' && k - 1 >= 0 && ss[k - 1] == '(')) {
					k++;
					j++;
					continue;
				}
				chain[j].residues.push_back(rna->chains[i1].residues[i2]);
				k++;
				i2++;
			}
		}
		rna->chains.clear();
		for (int i1 = 0; i1 < temp; i1++) {
			rna->chains.push_back(chain[i1]);
		}
		delete [] chain;
	}

	// mutate
	log("  mutate...");
	for (int i = 0; i < results->getLen(); i++) {
		RNA *rna = (*results)[i];
		mutate(rna);
	}

	// add P
	log("  add phosphate group...");
	for (int i = 0; i < results->getLen(); i++) {
		RNA *rna = results->at(i);
		addP(rna);
	}

	// assemble done
	ofile.open(logfile.c_str(), ios::app);
	ofile << Time::getTime() << "  done..." << endl << endl;
	ofile.close();
}

RNAs *Assemble::findRNA(loop *l, int counts) {
	static int iteration = 0;
	iteration++;

	// calculate number of loops to be created
	int loopCounts = l->getLoopCounts();
	double m_ = pow(counts, 1.0 / loopCounts);
	int m = int(ceil(m_));

	// find templates of itself
	RNAs *rnas = new RNAs;
	int a1 = l->s.getLen() - 1;
	int a2 = a1 + 1;
	if (l->s.head == NULL) {
		rnas = findLoop(l, m); 
	} else {
		RNA *rna = findHelix(&(l->s));
		if (l->head == NULL) {
			rnas->push(rna);
			return rnas;
		} else {
			RNAs *tempRNAs = findLoop(l, m);
			for (int i = 0; i < tempRNAs->getLen(); i++) {
				RNA *tempRNA = rna;
				tempRNA = new RNA(connect(*tempRNA, *(tempRNAs->at(i)), a1, a2));
				rnas->push(tempRNA); 
			}
			delete tempRNAs;
		}
		delete rna;
	}

//	int flag;
//	queue<int> a1_, a2_;
//	if (l->s.head == NULL) {
//		flag = -1;
//	} else {
//		flag = a1 - 2;
//	}
//	int n = 0;
//	for (res *tempRes = l->head; tempRes != NULL; tempRes = tempRes->next) {
//		if (tempRes->type != '&') {
//			flag++;
//		}
//		if (tempRes->type == '(') {
//			if (n % 2 == 1) {
//				if ((l->s.head == NULL && tempRes == l->head->next) || tempRes != l->head->next) {
//					a1_.push(flag);
//					a2_.push(flag + 1);
//				}
//			}
//			n++;
//		}
//	}

	// assemble
	int n = (l->s.head == NULL) ? 0 : (l->s.getLen() - mol->hinge_base_pair_num);
	int hinge_index = 0;
	for (loop *tempLoop = l->son; tempLoop != NULL; tempLoop = tempLoop->brother, hinge_index++) {
		loopCounts = tempLoop->getLoopCounts();
		RNAs *tempRNAs = findRNA(tempLoop, int(ceil(pow(double(m_), double(loopCounts)))));
		
		RNAs *tempRNAs2 = new RNAs;
//cerr << "(" << l->hinges[hinge_index].first << ", " << l->hinges[hinge_index].second << "), " << n << endl;
		for (int i = 0, aa = 0; i < rnas->getLen(); i++) {
			for (int j = 0; j < tempRNAs->getLen() && aa < counts; j++, aa++) {
				RNA *tempRNA2 = new RNA(connect(*((*rnas)[i]), *((*tempRNAs)[j]), 
				                        l->hinges[hinge_index].first + n,
				                        l->hinges[hinge_index].second + n));
				tempRNAs2->push(tempRNA2);
			}
		}

		n += tempRNAs->at(0)->getLen() - 4;

		delete rnas;
		delete tempRNAs;
		rnas = tempRNAs2;
		if (rnas->getLen() >= counts) {
			m_ = 1.0;
			m = 1;
		}
//		a1_.pop();
//		a2_.pop();
	}
	return rnas;
}

RNAs *Assemble::findLoop(loop *l, int counts) {
	if (l->head == NULL) {
		return NULL;
	}

//	// get loop info
//	string tempString;
//	for (res *r = l->head; r != NULL; r = r->next) {
//		if (r->type == '[' || r->type == ']') {
//			tempString += '.';
//		} else {
//			tempString += r->type;
//		}
//	}

	// find loop from loop library
	string l_seq = l->getSeq();
	string l_ss = l->getSS();
	string l_purified_ss;
	copy_if(l_ss.begin(), l_ss.end(), back_inserter(l_purified_ss), [](char c) {
		return c != '&';
	});
	string l_family = family;
	int hinge_num = l->hinges.size();

	stringstream strstr;
	string info_file;
	if (count_if(l_ss.begin(), l_ss.end(), [](char c) {
		return c == '[' || c == ']';
	})) {
		strstr << lib << "/info/pseudo_knot";
	} else {
		strstr << lib << "/info/loop-" << hinge_num;
	}
	strstr >> info_file;

//	int len = l->getLen();
//	int type = l->getType();
//	loopInfo *li = new loopInfo;
//	loopInfoList *liList = new loopInfoList;
//	string name = lib + "info/loop";
	ifstream ifile(info_file.c_str());
	if (!ifile) {
		cerr << "Assemble::findLoop error! Can't find '" << info_file << "'!" << endl;
		exit(1);
	}

	class LoopItem {
	public:
		string name;
		string seq;
		string ss;
		string family;
		int score;
	};
	LoopItem loop_item;

	vector<LoopItem> loop_sets;
	while (ifile >> loop_item.name >> loop_item.seq >> 
	       loop_item.ss >> loop_item.family) {
		if (loop_item.ss == l_ss) {
			loop_item.score = 0;

			// If this is a test case, then skip the templates from itself;
			// If it's not a test case, then increase the score by 2.
			if (loop_item.name.substr(5).compare(job.substr(0, 4)) == 0) {
				if (is_test) {
					continue;
				} else {
					loop_item.score += 5;
				}
			}

			// If this is a test case, then skip the templates from 
			// RNAs the family of which is the same;
			// If it's not a test case, then increase the score by 2.
			if (l_family == loop_item.family && l_family != "other") {
				if (is_test) {
					continue;
				} else {
					loop_item.score += 2;
				}
			}

			for (int i = 0; i < l_purified_ss.size(); i++) {
				if (l_seq[i] == loop_item.seq[i]) {
					if (set<char>{'(', ')', '[', ']'}.count(l_purified_ss[i])) {
						loop_item.score += 0.5;
					} else {
						loop_item.score += 1;
					}
				}
			}

			loop_sets.push_back(loop_item);
		}
	}
	ifile.close();
	sort(loop_sets.begin(), loop_sets.end(), [](
		const LoopItem &loop1, const LoopItem &loop2) {
			return loop1.score > loop2.score;
		}
	);

	ofstream ofile(logfile.c_str(), ios::app);
	ofile << "  for loop: ";
	ofile << l_seq << ' ' << l_ss << endl;

	RNAs *rnas = new RNAs;
	for (int i = 0; i < counts && i < loop_sets.size(); i++) {
		// get loop name
		string file_name = lib + "/loop/" + loop_sets[i].name + ".pdb";

		// log
		ofile << "    find: " << loop_sets[i].name << endl;
		ofile << "      sequence: " << loop_sets[i].seq << endl;
		ofile << "      secondary structure: " << loop_sets[i].ss << endl;
		ofile << "      family: " << loop_sets[i].family << endl;
		ofile << "      score: " << loop_sets[i].score << endl;

		RNA *rna = new RNA(file_name);
		rnas->push(rna);
	}

	if (loop_sets.size() < counts) {
		for (int i = 0; i < counts - loop_sets.size(); i++) {
			ofile << "    create loop..." << endl;
			ofile << "      method: Distance Geometry" << endl;
//			LoopModelling lm(l_seq, l_ss);
//			rnas->push(lm.run());
			LoopModelling2 lm;
			lm._view = _view;
			rnas->push(new RNA(lm(l_seq, l_ss)));
			ofile << "      final penalty energy: " << lm.dg.E << endl;
		}
	} 

	return rnas;
}

RNA *Assemble::findHelix(helix *s) {
	RNA *rna;
	bp *b;
	int i;
	string seq1, seq2;

	// if s is null, then return
	if (s->head == NULL) {
		return NULL;
	}

	// find helix info
	int len = s->getLen();
	string name = lib;
	name += "info/helix";
	ifstream ifile(name.c_str());
	helixInfo *hi = new helixInfo;
	helixInfoList *hiList = new helixInfoList;
	while (ifile >> hi->name >> hi->n >> hi->seq >> hi->ss >> hi->src) {
		if (len == hi->n) {
			int temp;
			double score;
			for (b = s->head, temp = 0, score = 0; b != NULL; b = b->next) {
				if (b->res1.name == hi->seq[temp]) {
					if (temp == 0 || temp == hi->n - 1) {
						score += 1 + 1 / (double) (2 * hi->n);
					} else {
						score += 1 / (double) (2 * hi->n);
					}
				}
				if (b->res2.name == hi->seq[2 * len - 1 - temp]) {
					if (temp == 0 || temp == hi->n - 1) {
						score += 1 + 1 / (double) (2 * hi->n);
					} else {
						score += 1 / (double) (2 * hi->n);
					}
				}
				temp++;
			}
			hiList->add(hi, score);
			hi = new helixInfo;
			if (hiList->getLen() == 5000) {
				break;
			}
			continue;
		}
		delete hi;
		hi = new helixInfo;
	}
	ifile.close();

	ofstream ofile(logfile.c_str(), ios::app);
	ofile << "  for helix: ";
	for (int i = 0; i < s->getLen(); i++) {
		ofile << '(';
	}
	for (int i = 0; i < s->getLen(); i++) {
		ofile << ')';
	}
	ofile << endl;

	if (hiList->head == NULL) {
		ofile << "    create helix..." << endl;
		ofile.close();
		seq1 = "";
		seq2 = "";
		for (b = s->head; b != NULL; b = b->next) {
			seq1 += b->res1.name;
			seq2 += b->res2.name;
		}
		for (i = seq2.size() - 1; i >= 0; i--) {
			seq1 += seq2[i];
		}
		return createHelix(seq1);
	} else {
		hi = hiList->head->hi;
		name = lib + "helix/";
		name += hi->name;
		name += ".pdb";
		ofile << "    find: " << name << endl;
		ofile << "      sequence: " << hi->seq << endl;
		ofile << "      score: " << hiList->head->score << endl;
		ofile.close();
		rna = new RNA(name);
		return rna;
	}
}

RNA *Assemble::createHelix(helix *s) {
	string seq1 = "";
	string seq2 = "";
	for (bp *b = s->head; b != NULL; b = b->next) {
		seq1 += b->res1.name;
		seq2 += b->res2.name;
	}
	for (int i = seq2.size() - 1; i >= 0; i--) {
		seq1 += seq2[i];
	}
	return createHelix(seq1);
}

RNA *Assemble::createHelix(string seq) {
	std::stringstream strstr;
	std::string file_name;

	if (seq.size() < 4) {
		std::cerr << "Assemble::createHelix error! The length of the helix"
		          << " to be create should not be less than 4!" << std::endl;
		exit(1);
	} else if (seq.size() % 2 == 1) {
		std::cerr << "Assemble::createHelix error! The length of the helix"
		          << " to be create should be even!" << std::endl;
		exit(1);
	} else if (seq.size() == 4 || seq.size() == 6) {
		strstr << lib << "/basepair/" << seq << ".pdb";
		strstr >> file_name;
		ifstream ifile(file_name.c_str());
		if (!ifile) {
			//strstr.str("");
			strstr.clear();
			strstr.seekp(0);
			strstr.seekg(0);
			strstr << lib << "/basepair/XXXXXX.pdb";
			strstr >> file_name;
		}
		ifile.close();
		return new RNA(file_name);
	} else {
		strstr << lib << "/basepair/" << seq.substr(0, 3) 
		       << seq.substr(seq.size() - 3, 3) << ".pdb";
		strstr >> file_name;
		ifstream ifile(file_name.c_str());
		if (!ifile) {
			//strstr.str("");
			strstr.clear();
			strstr.seekp(0);
			strstr.seekg(0);
			strstr << lib << "/basepair/XXXXXX.pdb";
			strstr >> file_name;
		}
		ifile.close();
		RNA *temp_rna_1 = new RNA(file_name);
		RNA *temp_rna_2 = createHelix(seq.substr(1, seq.size() - 2));
		RNA *return_rna = new RNA(connect(*temp_rna_1, *temp_rna_2, 2, 3));
		delete temp_rna_1;
		delete temp_rna_2;
		return return_rna;
	}
}

RNA *Assemble::createStrand(int n) {
	string name = lib;
	name += "longest_strand.pdb";
	RNA *oldRNA = new RNA(name.c_str());
	if (n > oldRNA->getLen()) {
		return NULL;
	}

	RNA *newRNA = new RNA;
	Chain *chain = new Chain;
	int num = 0;
	for (int i = 0; i < (int) oldRNA->chains.size() && num <= n; i++) {
		for (int j = 0; j < (int) oldRNA->chains[i].residues.size() && num <= n; j++) {
			num++;
			chain->push(oldRNA->chains[i].residues[j]);
		}
	}
	newRNA->push(chain);
	return newRNA;
}

void Assemble::addP(RNA *rna) {
	int i, j, k, temp;
	vector<Atom> va;
	Point *p = NULL, *q = NULL;
	Matr_ *a, *b, *c;
	double x1, x2, x, y1, y2, y, z1, z2, z;
	string name;
	Atom *atom;
	double r, r1, r2;

	for (i = 0; i < (int) rna->chains.size(); i++) {
		for (j = 0; j < (int) rna->chains[i].residues.size(); j++) {
			p = new Point[5];
			for (k = 0; k < (int) rna->chains[i].residues[j].atoms.size(); k++) {
				if (rna->chains[i].residues[j].atoms[k].name == "P") {
		      p[0].x = rna->chains[i].residues[j].atoms[k].x;
					p[0].y = rna->chains[i].residues[j].atoms[k].y;
			    p[0].z = rna->chains[i].residues[j].atoms[k].z;
			  } else if (rna->chains[i].residues[j].atoms[k].name == "O3*") {
		     	p[1].x = rna->chains[i].residues[j].atoms[k].x;
				 	p[1].y = rna->chains[i].residues[j].atoms[k].y;
					p[1].z = rna->chains[i].residues[j].atoms[k].z;
			  } else if (rna->chains[i].residues[j].atoms[k].name == "O5*") {
			    p[2].x = rna->chains[i].residues[j].atoms[k].x;
			    p[2].y = rna->chains[i].residues[j].atoms[k].y;
			    p[2].z = rna->chains[i].residues[j].atoms[k].z;
				} else if (rna->chains[i].residues[j].atoms[k].name == "C3*") {
			    p[3].x = rna->chains[i].residues[j].atoms[k].x;
		      p[3].y = rna->chains[i].residues[j].atoms[k].y;
					p[3].z = rna->chains[i].residues[j].atoms[k].z;
			  } else if (rna->chains[i].residues[j].atoms[k].name == "C5*") {
			    p[4].x = rna->chains[i].residues[j].atoms[k].x;
			    p[4].y = rna->chains[i].residues[j].atoms[k].y;
			    p[4].z = rna->chains[i].residues[j].atoms[k].z;
				}
			}
			for (k = 0, temp = 0; k < (int) rna->chains[i].residues[j].atoms.size(); k++) {
				name = rna->chains[i].residues[j].atoms[k].name;
				if (name == "P" || name == "O1P" || name == "O2P") {
					temp++;
				}
			}
			
			if (temp ==	3 && j != 0) {
				x1 = p[0].x - p[2].x; y1 = p[0].y - p[2].y; z1 = p[0].z - p[2].z;
				x2 = p[0].x - q[1].x; y2 = p[0].y - q[1].y; z2 = p[0].z - q[1].z;
	
				r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
			}

			if ((temp != 3 || (temp == 3 && (r1 > 2 || r2 > 2))) && j != 0) {
				atom = new Atom[3];
				atom[0].name = "P"; atom[1].name = "O1P"; atom[2].name = "O2P";
				x1 = q[1].x - q[3].x; y1 = q[1].y - q[3].y; z1 = q[1].z - q[3].z;
				x2 = p[2].x - p[4].x; y2 = p[2].y - p[4].y; z2 = p[2].z - p[4].z;
				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= r; y1 /= r; z1 /= r;
				r = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
				x2 /= r; y2 /= r; z2 /= r;
				x = (x1 + x2) / 2; y = (y1 + y2) / 2; z = (z1 + z2) / 2;
				x1 = p[2].x - q[1].x; y1 = p[2].y - q[1].y; z1 = p[2].z - q[1].z;
				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= r; y1 /= r; z1 /= r;
				r = -(x1 * x + y1 * y + z1 * z);
				x1 *= r; y1 *= r; z1 *= r;
				x += x1; y += y1; z += z1;
				r = sqrt(x * x + y * y + z * z);
				x /= r; y /= r; z /= r;
				atom[0].x = 0.884 * x + (q[1].x + p[2].x) / 2; atom[0].y = 0.884 * y + (q[1].y + p[2].y) / 2; atom[0].z = 0.884 * z + (q[1].z + p[2].z) / 2;
//				atom[0].x = (q[1].x + p[2].x) / 2; atom[0].y = (q[1].y + p[2].y) / 2; atom[0].z = (q[1].z + p[2].z) / 2;
				x1 = p[2].x - q[1].x; y1 = p[2].y - q[1].y; z1 = p[2].z - q[1].z;
				x2 = x; y2 = y; z2 = z;
				a = new Matr_(2, 2);
				b = new Matr_(2, 1);
				a->data[0][0] = x1; a->data[0][1] = y1; a->data[1][0] = x2; a->data[1][1] = y2;
				b->data[0][0] = -z1; b->data[1][0] = -z2;
				c = a->inverse()->multiply(b);
				x1 = c->data[0][0]; y1 = c->data[1][0]; z1 = 1;
				delete a;
				delete b;
				delete c;
				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= r; y1 /= r; z1 /= r;
				x1 *= 1.25; y1 *= 1.25; z1 *= 1.25;
				if ((p[4].x - p[2].x) * x1 + (p[4].y - p[2].y) * y1 + (p[4].z - p[2].z) * z1 > 0) {
					atom[1].x = atom[0].x + x1 + 0.884 * x;
					atom[1].y = atom[0].y + y1 + 0.884 * y;
					atom[1].z = atom[0].z + z1 + 0.884 * z;
					atom[2].x = atom[0].x - x1 + 0.884 * x;
					atom[2].y = atom[0].y - y1 + 0.884 * y;
					atom[2].z = atom[0].z - z1 + 0.884 * z;
				} else {
					atom[1].x = atom[0].x - x1 + 0.884 * x;
					atom[1].y = atom[0].y - y1 + 0.884 * y;
					atom[1].z = atom[0].z - z1 + 0.884 * z;
					atom[2].x = atom[0].x + x1 + 0.884 * x;
					atom[2].y = atom[0].y + y1 + 0.884 * y;
					atom[2].z = atom[0].z + z1 + 0.884 * z;
				}
				for (k = 0; k < (int) rna->chains[i].residues[j].atoms.size(); k++) {
					name = rna->chains[i].residues[j].atoms[k].name;
					if (name != "P" && name != "O1P" && name != "O2P") {
						va.push_back(rna->chains[i].residues[j].atoms[k]);
					}
				}
				rna->chains[i].residues[j].atoms.clear();
				rna->chains[i].residues[j].atoms.push_back(atom[0]);
				rna->chains[i].residues[j].atoms.push_back(atom[1]);
				rna->chains[i].residues[j].atoms.push_back(atom[2]);
				for (k = 0; k < (int) va.size(); k++) {
					rna->chains[i].residues[j].atoms.push_back(va[k]);
				}
			}
			if (q != NULL) {
				delete [] q;
			}
			q = p;
			va.clear();
		}
	}
}

void Assemble::mutate(RNA *r) {
	int i, j, k, k1, k2, k3, k_1, k_2, k_3, flag, temp;
	string resName1, resName2, baseName, name;
	double x, y, z;

	for (i = 0, flag = 0, temp = 0; i < (int) r->chains.size(); i++) {
		for (j = 0; j < (int) r->chains[i].residues.size(); j++, flag++) {
			if (flag >= (int) seq.size()) break;
			resName1 = r->chains[i].residues[j].name;
			if (seq[temp] != 'X') {
				resName2 = "";
				resName2 += seq[temp];
				if (r->chains[i].residues[j].name != resName2) {
					baseName = lib;
					baseName += "base/";
					baseName += resName2;
					baseName += ".pdb";
					RNA *rna = new RNA(baseName);

					// move rna1 to align the C1' atom 
					for (k = 0; k < (int) r->chains[i].residues[j].atoms.size(); k++) {
						if (r->chains[i].residues[j].atoms[k].name == "C1*") {
							k1 = k;
							x = r->chains[i].residues[j].atoms[k].x;
							y = r->chains[i].residues[j].atoms[k].y;
							z = r->chains[i].residues[j].atoms[k].z;
						}
						if (r->chains[i].residues[j].atoms[k].name == "C2") {
							k2 = k;
						}
					}
					for (k = 0; k < (int) r->chains[i].residues[j].atoms.size(); k++) {
						name = r->chains[i].residues[j].atoms[k].name;
						if ((name == "N9" && (resName1 == "A" || resName1 == "G")) || (name == "N1" && (resName1 == "U" || resName1 == "C"))) {
							k3 = k;
							goto out1;
						}
					}
					out1:
					r->move(-x, -y, -z);

					// move rna2 to align the C1' atom 
					for (k = 0; k < (int) rna->chains[0].residues[0].atoms.size(); k++) {
						if (rna->chains[0].residues[0].atoms[k].name == "C1*") {
							x = rna->chains[0].residues[0].atoms[k].x;
							y = rna->chains[0].residues[0].atoms[k].y;
							z = rna->chains[0].residues[0].atoms[k].z;
							k_1 = k;
						}
						if (rna->chains[0].residues[0].atoms[k].name == "C2") {
							k_2 = k;
						}
					}
					for (k = 0; k < (int) rna->chains[0].residues[0].atoms.size(); k++) {
						name = rna->chains[0].residues[0].atoms[k].name;
						if ((name == "N9" && (resName2 == "A" || resName2 == "G")) || (name == "N1" && (resName2 == "U" || resName2 == "C"))) {
							k_3 = k;
							goto out2;
						}
					}
					out2:
					rna->move(-x, -y, -z);

					// rotate rna1 with x-axis to make N atom on the x-z plane
					Matr_ *m = new Matr_(3, 3);
					x = r->chains[i].residues[j].atoms[k3].x;
					y = r->chains[i].residues[j].atoms[k3].y;
					z = r->chains[i].residues[j].atoms[k3].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						r->rotate(m);
					}
					delete m;
				
					// rotate rna2 with x-axis to make N atom on the x-z plane
					m = new Matr_(3, 3);
					x = rna->chains[0].residues[0].atoms[k_3].x;
					y = rna->chains[0].residues[0].atoms[k_3].y;
					z = rna->chains[0].residues[0].atoms[k_3].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						rna->rotate(m);
					}
					delete m;
				
					// rotate rna1 with y-axis to make N atom on the x-axis
					m = new Matr_(3, 3);
					x = r->chains[i].residues[j].atoms[k3].x;
					y = r->chains[i].residues[j].atoms[k3].y;
					z = r->chains[i].residues[j].atoms[k3].z;
					if (x * x + z * z != 0) {
						m->data[0][0] = x / sqrt(x * x + z * z);  m->data[0][1] = 0; m->data[0][2] = z / sqrt(x * x + z * z);
						m->data[1][0] = 0;                        m->data[1][1] = 1; m->data[1][2] = 0;
						m->data[2][0] = -z / sqrt(x * x + z * z); m->data[2][1] = 0; m->data[2][2] = x / sqrt(x * x + z * z);
						r->rotate(m);
					}
					delete m;

					// rotate rna2 with y-axis to make N atom on the x-axis
					m = new Matr_(3, 3);
					x = rna->chains[0].residues[0].atoms[k_3].x;
					y = rna->chains[0].residues[0].atoms[k_3].y;
					z = rna->chains[0].residues[0].atoms[k_3].z;
					if (x * x + z * z != 0) {
						m->data[0][0] = x / sqrt(x * x + z * z);  m->data[0][1] = 0; m->data[0][2] = z / sqrt(x * x + z * z);
						m->data[1][0] = 0;                        m->data[1][1] = 1; m->data[1][2] = 0;
						m->data[2][0] = -z / sqrt(x * x + z * z); m->data[2][1] = 0; m->data[2][2] = x / sqrt(x * x + z * z);
						rna->rotate(m);
					}
					delete m;

					// rotate rna1 with x-axis to make C2 atom on the x-z plane

					m = new Matr_(3, 3);
					x = r->chains[i].residues[j].atoms[k2].x;
					y = r->chains[i].residues[j].atoms[k2].y;
					z = r->chains[i].residues[j].atoms[k2].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						r->rotate(m);
					}
					delete m;

					// rotate rna2 with x-axis to make C2 atom on the x-z plane
					m = new Matr_(3, 3);
					x = rna->chains[0].residues[0].atoms[k_2].x;
					y = rna->chains[0].residues[0].atoms[k_2].y;
					z = rna->chains[0].residues[0].atoms[k_2].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						rna->rotate(m);
					}
					delete m;

					// change base
					vector<Atom> v;
					for (k = 0; k < (int) r->chains[i].residues[j].atoms.size(); k++) {
						name = r->chains[i].residues[j].atoms[k].name;
						if (name == "P" || name == "O1P" || name == "O2P" || name[name.size() - 1] == '*') {
							v.push_back(r->chains[i].residues[j].atoms[k]);
						}
					}
					r->chains[i].residues[j].atoms.clear();
					for (k = 0; k < (int) v.size(); k++) {
						r->chains[i].residues[j].atoms.push_back(v[k]);
					}
					for (k = 0; k < (int) rna->chains[0].residues[0].atoms.size(); k++) {
						if (rna->chains[0].residues[0].atoms[k].name != "C1*") {
							r->chains[i].residues[j].atoms.push_back(rna->chains[0].residues[0].atoms[k]);
						}
					}
					r->chains[i].residues[j].name = resName2;
					v.clear();
					delete rna;
				}
			}
			temp++;
		}
	}
}



