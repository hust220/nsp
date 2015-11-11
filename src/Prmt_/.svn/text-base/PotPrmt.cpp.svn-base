#include "PotPrmt.h"
#include "RNA.h"

PotPrmt::PotPrmt(char *str3) {
	string string3(str3);
	init(string3);
	theta = new double *[16];
	for (int i = 0; i < 16; i++) {
		theta[i] = new double[36];
	}
}

PotPrmt::PotPrmt(string str3) {
	init(str3);
}

void PotPrmt::init(string filename) {
	this->filename = filename;
	binWidth = 0.1;
	bins = int(1000 / binWidth);
}

void PotPrmt::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		cerr << str << endl;
		RNA rna(str);
		int len = rna.getLen();
		int type[len];
		Point **p = new Point *[len];
		for (int i = 0; i < len; i++) {
			p[i] = new Point[3];
		}
		int i = 0, j = 0, nj = 0, k = 0;
		while (i < rna.chains.size()) {
			if (k == 0) {
				if (rna[i][j].name == "A") {
					type[nj] = 0;
				} else if (rna[i][j].name == "U") {
					type[nj] = 1;
				} else if (rna[i][j].name == "G") {
					type[nj] = 2;
				} else if (rna[i][j].name == "C") {
					type[nj] = 3;
				}
			}
			Atom &atom = rna[i][j][k];
			if (type[nj] == 0) {
				if (atom.name == "C4") {
					p[nj][0].x = atom.x;
					p[nj][0].y = atom.y;
					p[nj][0].z = atom.z;
				} else if (atom.name == "N1") {
					p[nj][1].x = atom.x;
					p[nj][1].y = atom.y;
					p[nj][1].z = atom.z;
				} else if (atom.name == "N6") {
					p[nj][2].x = atom.x;
					p[nj][2].y = atom.y;
					p[nj][2].z = atom.z;
				}
			} else if (type[nj] == 1) {
				if (atom.name == "C6") {
					p[nj][0].x = atom.x;
					p[nj][0].y = atom.y;
					p[nj][0].z = atom.z;
				} else if (atom.name == "N3") {
					p[nj][1].x = atom.x;
					p[nj][1].y = atom.y;
					p[nj][1].z = atom.z;
				} else if (atom.name == "O4") {
					p[nj][2].x = atom.x;
					p[nj][2].y = atom.y;
					p[nj][2].z = atom.z;
				}
			} else if (type[nj] == 2) {
				if (atom.name == "C4") {
					p[nj][0].x = atom.x;
					p[nj][0].y = atom.y;
					p[nj][0].z = atom.z;
				} else if (atom.name == "N1") {
					p[nj][1].x = atom.x;
					p[nj][1].y = atom.y;
					p[nj][1].z = atom.z;
				} else if (atom.name == "O6") {
					p[nj][2].x = atom.x;
					p[nj][2].y = atom.y;
					p[nj][2].z = atom.z;
				}
			} else if (type[nj] == 3) {
				if (atom.name == "C6") {
					p[nj][0].x = atom.x;
					p[nj][0].y = atom.y;
					p[nj][0].z = atom.z;
				} else if (atom.name == "N3") {
					p[nj][1].x = atom.x;
					p[nj][1].y = atom.y;
					p[nj][1].z = atom.z;
				} else if (atom.name == "N4") {
					p[nj][2].x = atom.x;
					p[nj][2].y = atom.y;
					p[nj][2].z = atom.z;
				}
			}
			k++;
			if (rna.chains[i].residues[j].atoms.size() == k) {
				k = 0;
				j++;
				nj++;
				if (rna.chains[i].residues.size() == j) {
					j = 0;
					i++;
				}
			}
		}

		for (int i = 0; i < len; i++) {
			for (int j = i + 1; j < len; j++) {
				analyze(p[i], p[j], type[i], type[j], i, j);
			}
		}
		for (int i = 0; i < len; i++) {
			delete [] p[i];
		}
		delete [] p;
	}
	string temp[] = {"AA", "AU", "AG", "AC","UA", "UU", "UG", "UC","GA", "GU", "GG", "GC",  "CA", "CU", "CG", "CC"};
	string lib = getenv("RNA");
	for (int i = 0; i < 16; i++) {
		string f = lib + "test/" + temp[i];
		ofstream ofile(f.c_str());
		for (int j = 0; j < 36; j++) {
			ofile << j * 5 << '\t' << theta[i][j] << endl;
		}
		ofile.close();
	}
}

void PotPrmt::analyze(Point *nuc1, Point *nuc2, int type1, int type2, int ii, int jj) {
	Point *n1 = new Point[3];
	Point *n2 = new Point[3];
	for (int i = 0; i < 3; i++) {
		n1[i].x = nuc1[i].x;
		n1[i].y = nuc1[i].y;
		n1[i].z = nuc1[i].z;
		n2[i].x = nuc2[i].x;
		n2[i].y = nuc2[i].y;
		n2[i].z = nuc2[i].z;
	}
	Point center1(0.5 * (n1[0].x + n1[1].x), 0.5 * (n1[0].y + n1[1].y), 0.5 * (n1[0].z + n1[1].z));
	for (int i = 0; i < 3; i++) {
		n1[i].x -= center1.x;
		n1[i].y -= center1.y;
		n1[i].z -= center1.z;
		n2[i].x -= center1.x;
		n2[i].y -= center1.y;
		n2[i].z -= center1.z;
	}
	double r2, r, x_, y_, z_, c, s;
	r2 = n1[1].x * n1[1].x + n1[1].y * n1[1].y;
	if (r2 != 0) {
		r = sqrt(r2);
		c = n1[1].y / r;
		s = n1[1].x / r;
		for (int i = 0; i < 3; i++) {
			x_ = c * n1[i].x - s * n1[i].y;
			y_ = s * n1[i].x + c * n1[i].y;
			z_ = n1[i].z;
			n1[i].x = x_;
			n1[i].y = y_;
			n1[i].z = z_;
			x_ = c * n2[i].x - s * n2[i].y;
			y_ = s * n2[i].x + c * n2[i].y;
			z_ = n2[i].z;
			n2[i].x = x_;
			n2[i].y = y_;
			n2[i].z = z_;
		}
	}
	r2 = n1[1].y * n1[1].y + n1[1].z * n1[1].z;
	if (r2 != 0) {
		r = sqrt(r2);
		c = n1[1].y / r;
		s = -n1[1].z / r;
		for (int i = 0; i < 3; i++) {
			x_ = n1[i].x;
			y_ = c * n1[i].y - s * n1[i].z;
			z_ = s * n1[i].y + c * n1[i].z;
			n1[i].x = x_;
			n1[i].y = y_;
			n1[i].z = z_;
			x_ = n2[i].x;
			y_ = c * n2[i].y - s * n2[i].z;
			z_ = s * n2[i].y + c * n2[i].z;
			n2[i].x = x_;
			n2[i].y = y_;
			n2[i].z = z_;
		}
	}
	Point z(n1[2].y * n1[1].z - n1[2].z * n1[1].y, n1[2].z * n1[1].x - n1[2].x * n1[1].z, n1[2].x * n1[1].y - n1[2].y * n1[1].x);
	r2 = z.x * z.x + z.z * z.z;
	if (r2 != 0) {
		r = sqrt(r2);
		c = z.z / r;
		s = -z.x / r;
		for (int i = 0; i < 3; i++) {
			x_ = c * n1[i].x + s * n1[i].z;
			y_ = n1[i].y;
			z_ = -s * n1[i].x + c * n1[i].z;;
			n1[i].x = x_;
			n1[i].y = y_;
			n1[i].z = z_;
			x_ = c * n2[i].x + s * n2[i].z;
			y_ = n2[i].y;
			z_ = -s * n2[i].x + c * n2[i].z;;
			n2[i].x = x_;
			n2[i].y = y_;
			n2[i].z = z_;
		}
	}
	Point center2(0.5 * (n2[0].x + n2[1].x), 0.5 * (n2[0].y + n2[1].y), 0.5 * (n2[0].z + n2[1].z));
	double dZ2 = center2.z * center2.z;
	double dXY2 = center2.x * center2.x + center2.y * center2.y;
	Point a(n2[2].x - center2.x, n2[2].y - center2.y, n2[2].z - center2.z);
	Point b(n2[1].x - center2.x, n2[1].y - center2.y, n2[1].z - center2.z);
	Point z2(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
	r2 = z2.x * z2.x + z2.y * z2.y + z2.z * z2.z;
	double angle = 0;
	if (r2 != 0) {
		r = sqrt(r2);
		angle = acos(z2.z / r) / 3.1415927 * 180;
	}
	if (dZ2 <= 4  && dXY2 <= 49) {
		theta[type1 * 4 + type2][int(angle / 5.)]++;
	//	cout << "basepairing: " << ii + 1 << ':' << jj + 1 << ' ' << center2.x << ' ' << center2.y << endl;
	} else if (dZ2 >= 9 && dZ2 <= 25 && dXY2 <= 17) {
	//	cout << "stacking: " << ii + 1 << ':' << jj + 1 << ' ' << dZ2 << ' ' << dXY2 << ' ' << angle << endl;
	}
	return;
}






