#include "LoopModelling4.h"

LoopModelling4::LoopModelling4(string seq, string ss_) {
	this->seq = seq;
	int temp = 0;
	for (int i = 0; i < ss_.size(); i++) {
		if (ss_[i] == '.' || ss_[i] == '(' || ss_[i] == ')' || ss_[i] == '[' || ss_[i] == ']') {
			ss += ss_[i];
			temp++;
		} else if (ss_[i] == '&' && i > 0) {
			brk.push_back(temp - 1);
		}
	}
	int len1 = seq.size();
	int len2 = ss.size();
	if (len1 != len2) {
		cerr << "Sequence's length is not equal to second structure's length!" << endl;
		exit(1);
	}
	resLen = ss.size();
	len = ss.size() * 13;
	type = new int[resLen];
	for (int i = 0; i < resLen; i++) {
		if (seq[i] == 'A') {
			type[i] = 0;
		} else if (seq[i] == 'U') {
			type[i] = 1;
		} else if (seq[i] == 'G') {
			type[i] = 2;
		} else {
			type[i] = 3;
		}
	}
	bound = NULL;
	backbone = NULL;

	ct = NULL;
	ss2ct();
	chir = NULL;
	setChir();

	// nucleotide parameter
	string lib = getenv("RNA");
	string filename = lib + "prmt/LM4/nuc";
	ifstream ifile(filename.c_str());
	nucPrmt = new Matr_ *[4];
	for (int i = 0; i < 4; i++) {
		nucPrmt[i] = new Matr_(78, 2);
		for (int j = 0; j < 78; j++) {
			ifile >> nucPrmt[i]->data[j][0] >> nucPrmt[i]->data[j][1];
		}
	}
	ifile.close();
	// stacking parameter
	filename = lib + "prmt/LM4/stack";
	ifile.open(filename.c_str());
	stackPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		stackPrmt[i] = new Matr_(169, 2);
		for (int j = 0; j < 169; j++) {
			ifile >> stackPrmt[i]->data[j][0] >> stackPrmt[i]->data[j][1];
		}
	}
	ifile.close();
	// adjacent parameter
	filename = lib + "prmt/LM4/adj";
	ifile.open(filename.c_str());
	adjPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		adjPrmt[i] = new Matr_(169, 2);
		for (int j = 0; j < 169; j++) {
			ifile >> adjPrmt[i]->data[j][0] >> adjPrmt[i]->data[j][1];
		}
	}
	ifile.close();
	// basepair parameter
	filename = lib + "prmt/LM4/bp";
	ifile.open(filename.c_str());
	bpPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		bpPrmt[i] = new Matr_(169, 2);
		for (int j = 0; j < 169; j++) {
			ifile >> bpPrmt[i]->data[j][0] >> bpPrmt[i]->data[j][1];
		}
	}
	ifile.close();

	init();
	dg = new DG(bound);
}

void LoopModelling4::init() {
	delete bound;
	bound = new Matr_(len, len);
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (i == j) {
				bound->data[i][j] = 0;
			} else {
				bound->data[j][i] = 7;
				bound->data[i][j] = 999;
			}
		}
	}
	for (int i = 0; i < resLen; i++) {
		for (int j = 0; j < 13; j++) {
			for (int k = j + 1; k < 13; k++) {
				int n = j * 13 + k - (j + 1) * (j + 2) / 2;
				bound->data[i * 13 + j][i * 13 + k] = nucPrmt[type[i]]->data[n][1];
				bound->data[i * 13 + k][i * 13 + j] = nucPrmt[type[i]]->data[n][0];
			}
		}
	}
	for (int i = 0; i < resLen - 1; i++) {
		int type1 = type[i];
		int type2 = type[i + 1];
		for (int j = 0; j < 13; j++) {
			for (int k = 0; k < 13; k++) {
				int n = j * 13 + k;
				bound->data[i * 13 + j][(i + 1) * 13 + k] = adjPrmt[type1 * 4 + type2]->data[n][1];
				bound->data[(i + 1) * 13 + k][i * 13 + j] = adjPrmt[type1 * 4 + type2]->data[n][0];
			}
		}
	}
	if (ct == NULL) return;
	for (int i = 0; i < ct->row; i++) {
		int a = int(ct->data[i][0]);
		int b = int(ct->data[i][1]);
		int type1 = type[a];
		int type2 = type[b];
		if (a + 1 < resLen && ss[a + 1] != ')') {
			for (int j = 0; j < 13; j++) {
				for (int k = 0; k < 13; k++) {
					int n = j * 13 + k;
					bound->data[a * 13 + j][(a + 1) * 13 + k] = stackPrmt[type1 * 4 + type2]->data[n][1];
					bound->data[(a + 1) * 13 + k][a * 13 + j] = stackPrmt[type1 * 4 + type2]->data[n][0];
				}
			}
		}
		if (b + 1 < resLen) {
			for (int j = 0; j < 13; j++) {
				for (int k = 0; k < 13; k++) {
					int n = j * 13 + k;
					bound->data[b * 13 + j][(b + 1) * 13 + k] = stackPrmt[type1 * 4 + type2]->data[n][1];
					bound->data[(b + 1) * 13 + k][b * 13 + j] = stackPrmt[type1 * 4 + type2]->data[n][0];
				}
			}
		}
		for (int j = 0; j < 13; j++) {
			for (int k = 0; k < 13; k++) {
				int n = j * 13 + k;
				bound->data[a * 13 + j][b * 13 + k] = bpPrmt[type1 * 4 + type2]->data[n][1];
				bound->data[b * 13 + k][a * 13 + j] = bpPrmt[type1 * 4 + type2]->data[n][0];
			}
		}
	}
	for (int i = 0; i < (int) brk.size(); i++) {
		int a = brk.at(i);
		if (a + 1 < resLen) {
			for (int j = 0; j < 13; j++) {
				for (int k = 0; k < 13; k++) {
					int n = j * 13 + k;
					bound->data[a * 13 + j][(a + 1) * 13 + k] = 999;
					bound->data[(a + 1) * 13 + k][a * 13 + j] = 7;
				}
			}
		}
	}
}

LoopModelling4::~LoopModelling4() {
	delete bound;
	delete backbone;
	delete ct;
	delete dg;
	delete [] type;
	delete nucPrmt;
	delete bpPrmt;
	delete stackPrmt;
}

RNA *LoopModelling4::run() {
	delete backbone;
	backbone = dg->run();
	Point *p = NULL;
	RNA *rna = new RNA;
	Chain chain;
	rna->chains.push_back(chain);
	Point *o = NULL;
	for (int i = 0; i < resLen; i++) {
		delete [] p;
		p = new Point[13];
		for (int j = 0; j < 13; j++) {
			p[j].x = backbone->data[i * 13 + j][0];
			p[j].y = backbone->data[i * 13 + j][1];
			p[j].z = backbone->data[i * 13 + j][2];
		}
		double chi = (rand() % 1000) / 1000. * 360;
		Residue *res = buildNuc(p, chi, type[i]);
		rna->chains[0].residues.push_back(*res);
	}
	return rna;
}

void LoopModelling4::ss2ct() {
	int iTemp;
	vector<char> vcList;
	vector<int> viList;
	vector<char> vcList2;
	vector<int> viList2;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '(' || ss[i] == '[') {
			iTemp++;
		}
	}
	if (iTemp == 0) return;
	delete ct;
	ct = new Matr_(iTemp, 2);
	int flag = 0;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '(') {
			vcList.push_back('(');
			viList.push_back(i);
		} else if (ss[i] == '[') {
			vcList2.push_back('[');
			viList2.push_back(i);
		} else if (ss[i] == ')') {
			if (vcList.size() == 0) {
				cerr << "wrong secondary structure:\n" << ss << endl;
				exit(1);
			}
			ct->data[flag][0] = viList.back();
			ct->data[flag][1] = i;
			flag++;
			vcList.pop_back();
			viList.pop_back();
		} else if (ss[i] == ']') {
			if (vcList2.size() == 0) {
				cerr << "wrong secondary structure:\n" << ss << endl;
				exit(1);
			}
			ct->data[flag][0] = viList2.back();
			ct->data[flag][1] = i;
			flag++;
			vcList2.pop_back();
			viList2.pop_back();
		}
	}
	if (vcList.size() != 0 || vcList2.size() != 0) {
		cerr << "wrong secondary structure:\n" << ss << endl;
		exit(1);
	}
}

void LoopModelling4::setChir() {
	delete chir;
	chir = new Matr_(resLen, 5);
	for (int i = 0; i < resLen; i++) {
		chir->data[i][0] = i * 13 + 7;
		chir->data[i][1] = i * 13 + 9;
		chir->data[i][2] = i * 13 + 10;
		chir->data[i][3] = i * 13 + 11;
		chir->data[i][4] = -2.7;
	}
}

double LoopModelling4::at(Matr_ *bound, int i, int j, int len) {
    if (j >= len) {
        j -= len;
        return at(bound, j, i, len);
    }
    if (i >= len) {
        i -= len;
        return at(bound, j, i, len);
    }
    return bound->data[i][j];
}

void LoopModelling4::assign(Matr_ *bound, int i, int j, double d, int len) {
    if (j >= len) {
        j -= len;
        assign(bound, j, i, d, len);
        return;
    }
    if (i >= len) {
        i -= len;
        assign(bound, j, i, d, len);
        return;
    }
    bound->data[i][j] = d;
}

Residue *LoopModelling4::buildNuc(Point *p1, double chi, int type) {
	Point *p, *p3;
	int lent, len3;
	string f_name;
	ifstream ifile;
	string lib = getenv("RNA");
	if (type == 0) {
		len3 = 11;
		lent = 22;
		p = new Point[22];
		p3 = new Point[11];
		f_name = lib + "A";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 11; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 1) {
		len3 = 9;
		lent = 20;
		p = new Point[20];
		p3 = new Point[9];
		f_name = lib + "U";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 9; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 2) {
		len3 = 12;
		lent = 23;
		p = new Point[23];
		p3 = new Point[12];
		f_name = lib + "G";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 12; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 3) {
		len3 = 9;
		lent = 20;
		p = new Point[20];
		p3 = new Point[9];
		f_name = lib + "C";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 9; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else {
		cerr << "LoopModelling2::buildNuc error! The type must be one of 0, 1, 2, 3" << endl;
		exit(1);
	}

	// move p3 to origin
	double x_ = p3[0].x;
	double y_ = p3[0].y;
	double z_ = p3[0].z;
	for (int i = 0; i < len3; i++) {
		p3[i].x -= x_;
		p3[i].y -= y_;
		p3[i].z -= z_;
	}

	Point target(p1[12].x - p1[11].x, p1[12].y - p1[11].y, p1[12].z - p1[11].z);
	Point::coincide(p3, len3, p3[1], target);

	// normalVector
	Point *n1 = Point::normalVector(p1[6], p1[11], p1[9]);

	// dihedral
	double dih = Point::dihedral(n1, &(p3[0]), &(p3[1]), &(p3[2]));
	double delta_dih = chi - dih;

	// rotate
	Point origin;
	Point::rotate(p3, len3, origin, p3[1], delta_dih);
	for (int i = 0; i < len3; i++) {
		p3[i].x += p1[11].x;
		p3[i].y += p1[11].y;
		p3[i].z += p1[11].z;
	}

	for (int i = 0; i < 11; i++) {
		p[i].x = p1[i].x;
		p[i].y = p1[i].y;
		p[i].z = p1[i].z;
	}
	for (int i = 0; i < len3; i++) {
		p[i + 11].x = p3[i].x;
		p[i + 11].y = p3[i].y;
		p[i + 11].z = p3[i].z;
	}

	Residue *residue = new Residue(p, lent, type);
	return residue;
}




