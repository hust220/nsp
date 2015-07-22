#include "LM.h"

LM::LM(string seq, string ss_) {
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
	len = ss.size();
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

	init();
	dg = new DG(bound);
}

void LM::init() {
	delete bound;
	bound = new Matr_(len, len);
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (i == j) {
				bound->data[i][j] = 0;
			} else {
				bound->data[j][i] = 12;
				bound->data[i][j] = 999;
			}
		}
	}

	for (int i = 0; i < resLen - 1; i++) {
		bound->data[i][i + 1] = 11;
		bound->data[i + 1][i] = 13;
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
				bound->data[b * 5 + k][a * 5 + j] = 12;
			}
		}
	}

	if (ct != NULL) {
		for (int i = 0; i < ct->row; i++) {
			int a = int(ct->data[i][0]);
			int b = int(ct->data[i][1]);
			bound->data[a][b] = 15.6;
			bound->data[b][a] = 15.6;
			if (i + 1 < ct->row) {
				int c = int(ct->data[i + 1][0]);
				int d = int(ct->data[i + 1][1]);
				if (a + b == c + d && a - c == 1) {
					bound->data[c][a] = 6.3;
					bound->data[a][c] = 6.3;
					bound->data[b][d] = 6.3;
					bound->data[d][b] = 6.3;
					bound->data[c][b] = 17.5;
					bound->data[b][c] = 17.5;
					bound->data[a][d] = 12.5;
					bound->data[d][a] = 12.5;
				}
			}
			if (i + 2 < ct->row) {
				int c = int(ct->data[i + 2][0]);
				int d = int(ct->data[i + 2][1]);
				if (a + b == c + d && a - c == 2) {
					bound->data[c][a] = 12;
					bound->data[a][c] = 12;
					bound->data[b][d] = 12;
					bound->data[d][b] = 12;
					bound->data[c][b] = 19;
					bound->data[b][c] = 19;
					bound->data[a][d] = 11;
					bound->data[d][a] = 11;
				}
			}
		}
	}
}

RNA *LM::run() {
	delete coord;
	coord = dg->run();
	coord->print();
	Point *c = new Point[resLen];
	for (int i = 0; i < resLen; i++) {
		for (int j = 0; j < 3; j++) {
			c[i][j] = coord->data[i][j];
		}
	}
	NAST nast(c, resLen, ss);
	Point *c_ = nast.run();
	for (int i = 0; i < resLen; i++) {
		cout << c[i][0] << ' ' << c[i][1] << ' ' << c[i][2] << endl;
	}
//	mc(10000);
//	RNA *rna = newRNA();
	return NULL;
}

void LM::mc(int num) {
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

double LM::totEn() {
	double et = 0;
	int resLen = ss.size();
	for (int i = 0; i < resLen; i++) {
		et += energy(i, i, resLen - 1);
	}
	return et;
}

double LM::energy(int n, int left, int right) {
	double en = 0;
	for (int i = left; i <= right; i++) {
	}
	return en;
}

void LM::ss2ct() {
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

double LM::at(Matr_ *bound, int i, int j, int len) {
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

void LM::assign(Matr_ *bound, int i, int j, double d, int len) {
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

RNA *LM::newRNA() {
	return NULL;
	/*
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
	*/
}

Residue *LM::buildNuc(Matr_ *sixPt, int l, Point *base, int t) {
	return NULL;
	/*
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
	SupPos sp(m, n);
	SupPos sp_(m_, n);

	double rmsd1 = sp.getRMSD();
	double rmsd2 = sp_.getRMSD();
	Matr_ *rot;
	Point *c1;
	Point *c2;
	if (rmsd1 > rmsd2) {
		rot = sp_.getRot();
		c1 = sp_.getC1();
		c2 = sp_.getC2();
		delete cd;
		cd = cd_;
	} else {
		rot = sp.getRot();
		c1 = sp.getC1();
		c2 = sp.getC2();
	}
	delete m;
	delete m_;
	delete n;
	Matr_ *p1 = new Matr_(length - 1, 3);
	for (int i = 0; i < length - 1; i++) {
		for (int j = 0; j < 3; j++) {
			p1->data[i][j] = cd->data[i][j] - (*c1)[j];
		}
	}
	delete cd;
	delete c1;
	Matr_ *p2 = p1->multiply(rot);
	delete p1;
	delete rot;
	for (int i = 0; i < length - 1; i++) {
		for (int j = 0; j < 3; j++) {
			p2->data[i][j] += (*c2)[j];
		}
	}
	delete c2;
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
		cerr << "LM::buildNuc error!" << endl;
		exit(1);
	}
	Point *p = new Point[nucLen];
	for (int i = 0; i < length - 3; i++) {
		for (int j = 0; j < 3; j++) {
			p[i][j] = p2->data[i][j];
		}
	}
	delete p2;
	for (int i = length - 3; i < nucLen; i++) {
		for (int j = 0; j < 3; j++) {
			p[i][j] = base[i + 3 - length][j];
		}
	}
	Residue *res = new Residue(p, nucLen, t);
	delete [] p;
	return res;
	*/
}

Point **LM::cmpltBase() {
	return NULL;
	/*
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
	*/
}

