#include "LoopModelling.h"

LoopModelling::LoopModelling(string seq, string ss_) {
	this->seq = seq;
	int temp = 0;
	for (int i = 0; i < ss_.size(); i++) {
		if (ss_[i] == '.' || ss_[i] == '(' || ss_[i] == ')' || ss_[i] == '[' || ss_[i] == ']') {
			ss += ss_[i];
			temp++;
		} else if (ss_[i] == '&' && i > 0 && temp < seq.size()) {
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
	len = ss.size() * 5;
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
	coord = NULL;

	ct = NULL;
	ss2ct();

	// nucleotide parameter
	lib = getenv("RNA");
	string filename = lib + "prmt/LM/nuc";
	ifstream ifile(filename.c_str());
	nucPrmt = new Matr_ *[4];
	for (int i = 0; i < 4; i++) {
		nucPrmt[i] = new Matr_(10, 2);
		for (int j = 0; j < 10; j++) {
			ifile >> nucPrmt[i]->data[j][0] >> nucPrmt[i]->data[j][1];
		}
	}
	ifile.close();

	// adjacent parameter
	filename = lib + "prmt/LM/adj";
	ifile.open(filename.c_str());
	adjPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		adjPrmt[i] = new Matr_(25, 2);
		for (int j = 0; j < 25; j++) {
			ifile >> adjPrmt[i]->data[j][0] >> adjPrmt[i]->data[j][1];
		}
	}
	ifile.close();

	// stacking parameter
	filename = lib + "prmt/LM/stack";
	ifile.open(filename.c_str());
	stackPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		stackPrmt[i] = new Matr_(25, 2);
		for (int j = 0; j < 25; j++) {
			ifile >> stackPrmt[i]->data[j][0] >> stackPrmt[i]->data[j][1];
		}
	}
	ifile.close();

	// next parameter
	filename = lib + "prmt/LM/next";
	ifile.open(filename.c_str());
	nxtPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		nxtPrmt[i] = new Matr_(25, 2);
		for (int j = 0; j < 25; j++) {
			ifile >> nxtPrmt[i]->data[j][0] >> nxtPrmt[i]->data[j][1];
		}
	}
	ifile.close();

	// basepair parameter
	filename = lib + "prmt/LM/bp";
	ifile.open(filename.c_str());
	bpPrmt = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		bpPrmt[i] = new Matr_(25, 2);
		for (int j = 0; j < 25; j++) {
			ifile >> bpPrmt[i]->data[j][0] >> bpPrmt[i]->data[j][1];
		}
	}
	ifile.close();

	string nuc3 = lib + "prmt/LM/nuc3";
	ifile.open(nuc3.c_str());
	nucPrmt3 = new Matr_ *[4];
	lens = new int[4];
	lens[0] = 23;
	lens[1] = 21;
	lens[2] = 24;
	lens[3] = 21;
	for (int i = 0; i < 4; i++) {
		int length = lens[i] * (lens[i] - 1) / 2;
		nucPrmt3[i] = new Matr_(length, 2);
		for (int j = 0; j < length; j++) {
			ifile >> nucPrmt3[i]->data[j][0] >> nucPrmt3[i]->data[j][1];
		}
	}
	ifile.close();

	string nuc4 = lib + "prmt/LM/nuc4";
	ifile.open(nuc4.c_str());
	nucPrmt4 = new Matr_(91, 2);
	for (int i = 0; i < 91; i++) {
		ifile >> nucPrmt4->data[i][0] >> nucPrmt4->data[i][1];
	}
	ifile.close();

	init();
	dg = new DG(bound);
}

void LoopModelling::init() {
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
		for (int j = 0; j < 5; j++) {
			for (int k = j + 1; k < 5; k++) {
				int n = j * 5 + k - (j + 1) * (j + 2) / 2;
				bound->data[i * 5 + j][i * 5 + k] = nucPrmt[type[i]]->data[n][1];
				bound->data[i * 5 + k][i * 5 + j] = nucPrmt[type[i]]->data[n][0];
			}
		}
	}
	for (int i = 0; i < resLen - 1; i++) {
		int type1 = type[i];
		int type2 = type[i + 1];
		if (ss[i] == ss[i + 1] && ss[i] != '.') {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					int n = j * 5 + k;
					bound->data[i * 5 + j][(i + 1) * 5 + k] = stackPrmt[type1 * 4 + type2]->data[n][1];
					bound->data[(i + 1) * 5 + k][i * 5 + j] = stackPrmt[type1 * 4 + type2]->data[n][0];
				}
			}
		} else {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					int n = j * 5 + k;
					bound->data[i * 5 + j][(i + 1) * 5 + k] = adjPrmt[type1 * 4 + type2]->data[n][1];
					bound->data[(i + 1) * 5 + k][i * 5 + j] = adjPrmt[type1 * 4 + type2]->data[n][0];
				}
			}
		}
	}

	for (int i = 2; i < resLen; i++) {
		if (ss[i] != '.') {
			if (ss[i - 2] == ss[i - 1] && ss[i - 1] == ss[i]) {
				int type1 = type[i - 2];
				int type2 = type[i];
				for (int j = 0; j < 5; j++) {
					for (int k = 0; k < 5; k++) {
						int n = j * 5 + k;
						bound->data[(i - 2) * 5 + j][i * 5 + k] = nxtPrmt[type1 * 4 + type2]->data[n][1];
						bound->data[i * 5 + k][(i - 2) * 5 + j] = nxtPrmt[type1 * 4 + type2]->data[n][0];
					}
				}
			}
		}
	}

	for (int i = 0; i < (int) brk.size(); i++) {
		int a = brk[i];
		int b = a + 1;
		int type1 = type[a];
		int type2 = type[b];
		for (int j = 0; j < 5; j++) {
			for (int k = 0; k < 5; k++) {
				int n = j * 5 + k;
				bound->data[a * 5 + j][b * 5 + k] = 999;
				bound->data[b * 5 + k][a * 5 + j] = 7;
			}
		}
	}

	if (ct == NULL) return;
	for (int i = 0; i < ct->row; i++) {
		int a = int(ct->data[i][0]);
		int b = int(ct->data[i][1]);
		int type1 = type[a];
		int type2 = type[b];
		for (int j = 0; j < 5; j++) {
			for (int k = 0; k < 5; k++) {
				int n = j * 5 + k;
				bound->data[a * 5 + j][b * 5 + k] = bpPrmt[type1 * 4 + type2]->data[n][1];
				bound->data[b * 5 + k][a * 5 + j] = bpPrmt[type1 * 4 + type2]->data[n][0];
			}
		}
	}
}

RNA *LoopModelling::run() {
	delete coord;
	coord = dg->run();
//	coord->print();
//	mc(10000);
	RNA *rna = newRNA();
	return rna;
}

void LoopModelling::mc(int num) {
	MRand mr;
	int resLen = ss.size();
	double et = totEn();
	Matr_ *mMin = new Matr_(coord);
	double dMin = et;
	for (int i = 0; i < num; i++) {
		if (i % int(num / 100.) == 0) {
			cerr << i << ' ' << et << ' ' << dMin << endl;
		}
		int n = int(mr.run() * len);
		double eo = energy(n / 5, 0, resLen - 1);
		double x_ = coord->data[n][0];
		double y_ = coord->data[n][1];
		double z_ = coord->data[n][2];
		coord->data[n][0] += (mr.run() - 0.5) * 2;
		coord->data[n][1] += (mr.run() - 0.5) * 2;
		coord->data[n][2] += (mr.run() - 0.5) * 2;
		double en = energy(n / 5, 0, resLen - 1);
		if (fabs(eo) > 1.e-9 && mr.run() > exp(-(en - eo) / eo)) {
			coord->data[n][0] = x_;
			coord->data[n][0] = x_;
			coord->data[n][0] = x_;
		} else {
			et += en - eo;
			if (dMin > et) {
				for (int ii = 0; ii < mMin->row; ii++) {
					for (int jj = 0; jj < mMin->col; jj++) {
						mMin->data[ii][jj] = coord->data[ii][jj];
					}
				}
			}
		}
	}
	for (int ii = 0; ii < mMin->row; ii++) {
		for (int jj = 0; jj < mMin->col; jj++) {
			coord->data[ii][jj] = mMin->data[ii][jj];
		}
	}
}

double LoopModelling::totEn() {
	double et = 0;
	int resLen = ss.size();
	for (int i = 0; i < resLen; i++) {
		et += energy(i, i, resLen - 1);
	}
	return et;
}

double LoopModelling::energy(int n, int left, int right) {
	double en = 0;
	for (int i = left; i <= right; i++) {
	}
	return en;
}

void LoopModelling::ss2ct() {
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

double LoopModelling::at(Matr_ *bound, int i, int j, int len) {
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

void LoopModelling::assign(Matr_ *bound, int i, int j, double d, int len) {
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

RNA *LoopModelling::newRNA() {
	RNA *rna = new RNA;
	Chain *chain = new Chain;
	Point **base = cmpltBase();
	for (int i = 0; i < resLen; i++) {
		int l = 4;
		if (i == resLen - 1 || (ss[i] == '(' && ss[i + 1] == ')')) {
			l = 3;
		}
		Matr_ *sixPt = new Matr_(l, 3);
		for (int k = 0; k < 3; k++) {
			sixPt->data[0][k] = coord->data[i * 5][k];
			sixPt->data[1][k] = base[i][0][k];
			sixPt->data[2][k] = base[i][1][k];
			if (l == 5) {
				sixPt->data[3][k] = coord->data[i * 5 + 5][k];
			}
		}
		Residue *res = buildNuc(sixPt, l, base[i], type[i]);
		chain->push(res);
		int flag = 0;
		for (int j = 0; j < (int) brk.size(); j++) {
			if (brk[j] == i) {
				flag = 1;
				break;
			}
		}
		if ((ss[i] == '(' && ss[i + 1] == ')') || flag == 1) {
			rna->push(*chain);
			delete chain;
			chain = new Chain;
		}
		delete res;
		delete sixPt;
	}
	rna->push(*chain);
	delete chain;
	return rna;
}

Residue *LoopModelling::buildNuc(Matr_ *sixPt, int l, Point *base, int t) {
	int length = 14;
	Matr_ *b = new Matr_(length, length);
	for (int i = 0; i < length; i++) {
		for (int j = i; j < length; j++) {
			if (i == j) {
				b->data[i][j] = 0;
				continue;
			}
			b->data[i][j] = nucPrmt4->data[i * length + j - (i + 1) * (i + 2) / 2][1];
			b->data[j][i] = nucPrmt4->data[i * length + j - (i + 1) * (i + 2) / 2][0];
		}
	}

	int flag[4];
	flag[0] = 0;
	flag[1] = 11;
	flag[2] = 12;
	flag[3] = 13;
	for (int i = 0; i < l; i++) {
		for (int j = i + 1; j < l; j++) {
			Point p1(sixPt->data[i][0], sixPt->data[i][1], sixPt->data[i][2]);
			Point p2(sixPt->data[j][0], sixPt->data[j][1], sixPt->data[j][2]);
			double distance = p1.dist(p2);
			double l = b->data[flag[i]][flag[j]];
			double u = b->data[flag[j]][flag[i]];
			/*
			if (distance > u) {
				distance = u;
			} else if (distance < l) {
				distance = l;
			}
			*/
			b->data[flag[i]][flag[j]] = distance;
			b->data[flag[j]][flag[i]] = distance;
		}
	}

	DG dg_(b);
	Matr_ *cd = dg_.run();
	Matr_ *cd_ = new Matr_(cd);
	for (int i = 0; i < length; i++) {
		cd_->data[i][1] = -cd_->data[i][1];
	}
	delete b;
	Matr_ *m = new Matr_(l, 3);
	Matr_ *m_ = new Matr_(l, 3);
	Matr_ *n = new Matr_(l, 3);
	for (int i = 0; i < l; i++) {
		for (int j = 0; j < 3; j++) {
			m->data[i][j] = cd->data[flag[i]][j];
			m_->data[i][j] = cd_->data[flag[i]][j];
			n->data[i][j] = sixPt->data[i][j];
		}
	}
	SupPos sp, sp_;
	sp(m, n);
	sp_(m_, n);

	double rmsd1 = sp.rmsd;
	double rmsd2 = sp_.rmsd;
	Matrix3f rot;
	Point c1;
	Point c2;
	if (rmsd1 > rmsd2) {
		rot = sp_.rot;
		c1 = sp_.c1;
		c2 = sp_.c2;
		delete cd;
		cd = cd_;
	} else {
		rot = sp.rot;
		c1 = sp.c1;
		c2 = sp.c2;
	}
	delete m;
	delete m_;
	delete n;
	MatrixXf p1(length - 1, 3);
	for (int i = 0; i < length - 1; i++) {
		for (int j = 0; j < 3; j++) {
			p1(i, j) = cd->data[i][j] - c1[j];
		}
	}
	delete cd;
	MatrixXf p2 = p1 * rot;
	for (int i = 0; i < length - 1; i++) {
		for (int j = 0; j < 3; j++) {
			p2(i, j) += c2[j];
		}
	}
	int nucLen;
	if (t == 0) {
		nucLen = 22;
	} else if (t == 1) {
		nucLen = 20;
	} else if (t == 2) {
		nucLen = 23;
	} else if (t == 3) {
		nucLen = 20;
	} else {
		cerr << "LoopModelling::buildNuc error!" << endl;
		exit(1);
	}
	Point *p = new Point[nucLen];
	for (int i = 0; i < length - 3; i++) {
		for (int j = 0; j < 3; j++) {
			p[i][j] = p2(i, j);
		}
	}
	for (int i = length - 3; i < nucLen; i++) {
		for (int j = 0; j < 3; j++) {
			p[i][j] = base[i + 3 - length][j];
		}
	}
	Residue *res = new Residue(p, nucLen, t);
	delete [] p;
	return res;
/*
	double x[3];
	x[0] = cd->data[0][0];
	x[1] = cd->data[0][1];
	x[2] = cd->data[0][2];
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < 3; j++) {
			cd->data[i][j] -= x[j];
		}
	}

	Point *p = new Point[length];
	for (int i = 0; i < length; i++) {
		p[i].x = cd->data[i][0];
		p[i].y = cd->data[i][1];
		p[i].z = cd->data[i][2];
	}
	Point target(sixPt->data[1][0] - sixPt->data[0][0], 
				 sixPt->data[1][1] - sixPt->data[0][1],
				 sixPt->data[1][2] - sixPt->data[0][2]);
	Point::coincide(p, length, p[flag[1]], target);

	Point p1(p[flag[2]].x, p[flag[2]].y, p[flag[2]].z);
	if (l == 6) {
		p1.x = p[flag[5]].x;
		p1.y = p[flag[5]].y;
		p1.z = p[flag[5]].z;
	}
	Point p2(0, 0, 0);
	Point p3(p[flag[1]].x, p[flag[1]].y, p[flag[1]].z);
	Point p4(sixPt->data[2][0] - sixPt->data[0][0],
	         sixPt->data[2][1] - sixPt->data[0][1],
			 sixPt->data[2][2] - sixPt->data[0][2]);
	if (l == 6) {
		p4.x = sixPt->data[5][0] - sixPt->data[0][0];
		p4.y = sixPt->data[5][0] - sixPt->data[0][0];
		p4.z = sixPt->data[5][0] - sixPt->data[0][0];
	}
	double ang = Point::dihedral(p1, p2, p3, p4);
	Point origin;
	Point::rotate(p, length, origin, p[flag[1]], ang);
	for (int i = 0; i < length; i++) {
		p[i].x += sixPt->data[0][0];
		p[i].y += sixPt->data[0][1];
		p[i].z += sixPt->data[0][2];
	}
	Residue *res = new Residue(p, length - 1, t);

	return res;
*/
}

Point **LoopModelling::cmpltBase() {
	Point **base = new Point *[resLen];
	ifstream ifile;
	string filename;

	Point *bsA = new Point[11];
	filename = lib + "base/A";
	ifile.open(filename.c_str());
	for (int i = 0; i < 11; i++) {
		ifile >> bsA[i][0] >> bsA[i][1] >> bsA[i][2];
	}
	ifile.close();
	
	Point *bsU = new Point[9];
	filename = lib + "base/U";
	ifile.open(filename.c_str());
	for (int i = 0; i < 9; i++) {
		ifile >> bsU[i][0] >> bsU[i][1] >> bsU[i][2];
	}
	ifile.close();
	
	Point *bsG = new Point[12];
	filename = lib + "base/G";
	ifile.open(filename.c_str());
	for (int i = 0; i < 12; i++) {
		ifile >> bsG[i][0] >> bsG[i][1] >> bsG[i][2];
	}
	ifile.close();
	
	Point *bsC = new Point[9];
	filename = lib + "base/C";
	ifile.open(filename.c_str());
	for (int i = 0; i < 9; i++) {
		ifile >> bsC[i][0] >> bsC[i][1] >> bsC[i][2];
	}
	ifile.close();


	for (int i = 0; i < resLen; i++) {
		Point *temp;
		int flag[3];
		int nucLen;
		if (type[i] == 0) {
			nucLen = 11;
			temp = new Point[11];
			for (int j = 0; j < 11; j++) {
				for (int k = 0; k < 3; k++) {
					temp[j][k] = bsA[j][k];
				}
			}
			flag[0] = 10;
			flag[1] = 7;
			flag[2] = 6;
		} else if (type[i] == 1) {
			nucLen = 9;
			temp = new Point[9];
			for (int j = 0; j < 9; j++) {
				for (int k = 0; k < 3; k++) {
					temp[j][k] = bsU[j][k];
				}
			}
			flag[0] = 8;
			flag[1] = 4;
			flag[2] = 6;
		} else if (type[i] == 2) {
			nucLen = 12;
			temp = new Point[12];
			for (int j = 0; j < 12; j++) {
				for (int k = 0; k < 3; k++) {
					temp[j][k] = bsG[j][k];
				}
			}
			flag[0] = 11;
			flag[1] = 7;
			flag[2] = 6;
		} else if (type[i] == 3) {
			nucLen = 9;
			temp = new Point[9];
			for (int j = 0; j < 9; j++) {
				for (int k = 0; k < 3; k++) {
					temp[j][k] = bsC[j][k];
				}
			}
			flag[0] = 8;
			flag[1] = 4;
			flag[2] = 6;
		}
		double x_ = temp[flag[0]][0];
		double y_ = temp[flag[0]][1];
		double z_ = temp[flag[0]][2];
		for (int j = 0; j < nucLen; j++) {
			temp[j][0] -= x_;
			temp[j][1] -= y_;
			temp[j][2] -= z_;
		}
		Point o_(temp[flag[1]][0], temp[flag[1]][1], temp[flag[1]][2]);
		Point n_(coord->data[i * 5 + 3][0] - coord->data[i * 5 + 2][0], 
		         coord->data[i * 5 + 3][1] - coord->data[i * 5 + 2][1],
				 coord->data[i * 5 + 3][2] - coord->data[i * 5 + 2][2]);
		Point::coincide(temp, nucLen, o_, n_);
		Point p1_(temp[flag[2]][0], temp[flag[2]][1], temp[flag[2]][2]);
		Point p2_(0, 0, 0);
		Point p3_(temp[flag[1]][0], temp[flag[1]][1], temp[flag[1]][2]);
		Point p4_(coord->data[i * 5 + 4][0] - coord->data[i * 5 + 2][0], 
		          coord->data[i * 5 + 4][1] - coord->data[i * 5 + 2][1],
				  coord->data[i * 5 + 4][2] - coord->data[i * 5 + 2][2]);
		double angle = Point::dihedral(p1_, p2_, p3_, p4_);
		Point origin_;
		Point::rotate(temp, nucLen, p2_, p3_, angle);
		x_ = coord->data[i * 5 + 2][0];
		y_ = coord->data[i * 5 + 2][1];
		z_ = coord->data[i * 5 + 2][2];
		for (int j = 0; j < nucLen; j++) {
			temp[j][0] += x_;
			temp[j][1] += y_;
			temp[j][2] += z_;
		}
		base[i] = temp;
	}
	delete [] bsA;
	delete [] bsU;
	delete [] bsG;
	delete [] bsC;
	return base;
}

