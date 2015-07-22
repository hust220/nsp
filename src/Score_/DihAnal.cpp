#include "DihAnal.h"

namespace jian {

DihAnal::DihAnal(double n) {
	p = NULL;
	o5_ = NULL;
	c5_ = NULL;
	c4_ = NULL;
	c3_ = NULL;
	o3_ = NULL;
	c2_ = NULL;
	c1_ = NULL;
	b1 = NULL;
	b2 = NULL;
	obsProb = NULL;
	refProb = NULL;
	len = 0;
	score = 0;
	interval = n;
	bins = int(360 / interval);
	obsParm = new int[7 * bins];
	for (int i = 0; i < 7 * bins; i++) {
		obsParm[i] = 0;
	}
	refParm = new int[bins];
	for (int i = 0; i < bins; i++) {
		refParm[i] = 0;
	}
}

DihAnal::~DihAnal() {
	delPoints();
	delete [] obsParm;
	delete [] obsProb;
	delete [] refParm;
	delete [] refProb;
}

void DihAnal::delPoints() {
	for (int i = 0; i < len; i++) {
		delete p[i];
		delete o5_[i];
		delete c5_[i];
		delete c4_[i];
		delete c3_[i];
		delete o3_[i];
		delete c2_[i];
		delete c1_[i];
		delete b1[i];
		delete b2[i];
	}
	delete [] p;
	delete [] o5_;
	delete [] c5_;
	delete [] c4_;
	delete [] c3_;
	delete [] o3_;
	delete [] c2_;
	delete [] c1_;
	delete [] b1;
	delete [] b2;
}

void DihAnal::initPoints(int len) {
	if (p != NULL) {
		delPoints();
	}

	p = new Point *[len];
	o5_ = new Point *[len];
	c5_ = new Point *[len];
	c4_ = new Point *[len];
	c3_ = new Point *[len];
	o3_ = new Point *[len];
	c2_ = new Point *[len];
	c1_ = new Point *[len];
	b1 = new Point *[len];
	b2 = new Point *[len];
	for (int i = 0; i < len; i++) {
		p[i] = NULL;
		o5_[i] = NULL;
		c5_[i] = NULL;
		c4_[i] = NULL;
		c3_[i] = NULL;
		o3_[i] = NULL;
		c2_[i] = NULL;
		c1_[i] = NULL;
		b1[i] = NULL;
		b2[i] = NULL;
	}
}

void DihAnal::readRNA(const RNA &r) {
	int length = r.res_nums();
	initPoints(length);
	len = length;

	for (int n = 0, i = 0; i < (int) r.chains.size(); i++) {
		for (int j = 0; j < (int) r.chains[i].residues.size(); j++, n++) {
			string resName = r.chains[i].residues[j].name;
			for (int k = 0; k < (int) r.chains[i].residues[j].atoms.size(); k++) {
				string name = r.chains[i].residues[j].atoms[k].name;
				double x, y, z;
				if (name == "P") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					p[n] = new Point(x, y, z);
				} else if (name == "O5*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					o5_[n] = new Point(x, y, z);
				} else if (name == "C5*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					c5_[n] = new Point(x, y, z);
				} else if (name == "C4*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					c4_[n] = new Point(x, y, z);
				} else if (name == "C3*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					c3_[n] = new Point(x, y, z);
				} else if (name == "O3*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					o3_[n] = new Point(x, y, z);
				} else if (name == "C2*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					c2_[n] = new Point(x, y, z);
				} else if (name == "C1*") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					c1_[n] = new Point(x, y, z);
				} else if ((resName == "A" || resName == "G") && name == "N9") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					b1[n] = new Point(x, y, z);
				} else if ((resName == "A" || resName == "G") && name == "C8") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					b2[n] = new Point(x, y, z);
				} else if ((resName == "U" || resName == "C") && name == "N1") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					b1[n] = new Point(x, y, z);
				} else if ((resName == "U" || resName == "C") && name == "C6") {
					x = r.chains[i].residues[j].atoms[k].x;
					y = r.chains[i].residues[j].atoms[k].y;
					z = r.chains[i].residues[j].atoms[k].z;
					b2[n] = new Point(x, y, z);
				}
			}
		}
	}
}

void DihAnal::train() {
	double temp;
	for (int i = 0; i < len; i++) {
		/* alpha */
		if (i >= 1 && p[i] != NULL && o3_[i - 1]->dist(o5_[i]) < 3) {
			temp = Point::dihedral(o3_[i - 1], p[i], o5_[i], c5_[i]);
			obsParm[int(temp / interval)]++;
			refParm[int(temp / interval)]++;
		}
		/* beta */
		if (p[i] != NULL) {
			temp = Point::dihedral(p[i], o5_[i], c5_[i], c4_[i]);
			obsParm[bins + int(temp / interval)]++;
			refParm[int(temp / interval)]++;
		}
		/* gamma */
		temp = Point::dihedral(o5_[i], c5_[i], c4_[i], c3_[i]);
		obsParm[bins * 2 + int(temp / interval)]++;
		refParm[int(temp / interval)]++;
		/* delta */
		temp = Point::dihedral(c5_[i], c4_[i], c3_[i], o3_[i]);
		obsParm[bins * 3 + int(temp / interval)]++;
		refParm[int(temp / interval)]++;
		/* epsilon */
		if (i + 1 < len && p[i + 1] != NULL && o3_[i]->dist(o5_[i + 1]) < 3) {
			temp = Point::dihedral(c4_[i], c3_[i], o3_[i], p[i + 1]);
			obsParm[bins * 4 + int(temp / interval)]++;
			refParm[int(temp / interval)]++;
		}
		/* zeta */
		if (i + 1 < len && p[i + 1] != NULL && o3_[i]->dist(o5_[i + 1]) < 3) {
			temp = Point::dihedral(c3_[i], o3_[i], p[i + 1], o5_[i + 1]);
			obsParm[bins * 5 + int(temp / interval)]++;
			refParm[int(temp / interval)]++;
		}
		/* chi */
		temp = Point::dihedral(c2_[i], c1_[i], b1[i], b2[i]);
		obsParm[bins * 6 + int(temp / interval)]++;
		refParm[int(temp / interval)]++;
	}
}

void DihAnal::readParm(char *filename) {
	string str(filename);
	readParm(str);
}

void DihAnal::readParm(string filename) {
	ifstream ifile(filename.c_str());
	int temp;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < bins; j++) {
			ifile >> temp;
			obsParm[i * bins + j] += temp;
		}
	}
	for (int i = 0; i < bins; i++) {
		ifile >> temp;
		refParm[i] += temp;
	}
	ifile.close();
	initProb();
}

double DihAnal::scoring() {
	score = 0;
	int n = 0;
	double temp, a, b;
	for (int i = 0; i < len; i++) {
		/* alpha */
		if (i >= 1 && p[i] != NULL && o3_[i - 1]->dist(o5_[i]) < 3) {
			temp = Point::dihedral(o3_[i - 1], p[i], o5_[i], c5_[i]);
			a = obsProb[int(temp / interval)];
			if (a == 0) {
				score -= 0;
			} else {
				b = refProb[int(temp / interval)];
				score -= log(a / b);
			}
			n++;
		}
		/* beta */
		if (p[i] != NULL) {
			temp = Point::dihedral(p[i], o5_[i], c5_[i], c4_[i]);
			a = obsProb[bins + int(temp / interval)];
			if (a == 0) {
				score -= 0;
			} else {
				b = refProb[int(temp / interval)];
				score -= log(a / b);
			}
			n++;
		}
		/* gamma */
		temp = Point::dihedral(o5_[i], c5_[i], c4_[i], c3_[i]);
		a = obsProb[bins * 2 + int(temp / interval)];
		if (a == 0) {
			score -= 0;
		} else {
			b = refProb[int(temp / interval)];
			score -= log(a / b);
		}
		n++;
		/* delta */
		temp = Point::dihedral(c5_[i], c4_[i], c3_[i], o3_[i]);
		a = obsProb[bins * 3 + int(temp / interval)];
		if (a == 0) {
			score -= 0;
		} else {
			b = refProb[int(temp / interval)];
			score -= log(a / b);
		}
		n++;
		/* epsilon */
		if (i + 1 < len && p[i + 1] != NULL && o3_[i]->dist(o5_[i + 1]) < 3) {
			temp = Point::dihedral(c4_[i], c3_[i], o3_[i], p[i + 1]);
			a = obsProb[bins * 4 + int(temp / interval)];
			if (a == 0) {
				score -= 0;
			} else {
				b = refProb[int(temp / interval)];
				score -= log(a / b);
			}
			n++;
		}
		/* zeta */
		if (i + 1 < len && p[i + 1] != NULL && o3_[i]->dist(o5_[i + 1]) < 3) {
			temp = Point::dihedral(c3_[i], o3_[i], p[i + 1], o5_[i + 1]);
			a = obsProb[bins * 5 + int(temp / interval)];
			if (a == 0) {
				score -= 0;
			} else {
				b = refProb[int(temp / interval)];
				score -= log(a / b);
			}
			n++;
		}
		/* chi */
		temp = Point::dihedral(c2_[i], c1_[i], b1[i], b2[i]);
		a = obsProb[bins * 6 + int(temp / interval)];
		if (a == 0) {
			score -= 0;
		} else {
			b = refProb[int(temp / interval)];
			score -= log(a / b);
		}
		n++;
	}
	score = score / n;
	return score;
}

void DihAnal::initProb() {
	if (obsProb == NULL) {
		obsProb = new double[7 * bins];
	}
	if (refProb == NULL) {
		refProb = new double[bins];
	}

	int *temp = new int[7];
	for (int i = 0; i < 7; i++) {
		temp[i] = 0;
	}
	int total = 0;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < bins; j++) {
			temp[i] += obsParm[i * bins + j];
		}
		total += temp[i];
		for (int j = 0; j < bins; j++) {
			obsProb[i * bins + j] = double(obsParm[i * bins + j]) / temp[i];
		}
	}
	for (int i = 0; i < bins; i++) {
		refProb[i] = double(refParm[i]) / total;
	}
	delete [] temp;
}

void DihAnal::printParm() {
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < bins; j++) {
			cout << obsParm[i * bins + j] << ' ';
		}
		cout << endl;
	}
	for (int i = 0; i < bins; i++) {
		cout << refParm[i] << ' ';
	}
	cout << endl;
}

void DihAnal::printProb() {
	if (obsProb == NULL) {
		obsProb = new double[7 * bins];
	}
	if (refProb == NULL) {
		refProb = new double[bins];
	}

	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < bins; j++) {
			cout << obsProb[i * bins + j] << ' ';
		}
		cout << endl;
	}
	for (int i = 0; i < bins; i++) {
		cout << refProb[i] << ' ';
	}
	cout << endl;
}

double DihAnal::getScore() {
	return score;
}

} /// namespace jian

