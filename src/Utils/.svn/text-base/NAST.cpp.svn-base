#include "NAST.h"

NAST::NAST(Point *c, int l, string s, Matr_ *t) {
	len = l;
	coord = new Point[len];
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < 3; j++) {
			coord[i][j] = c[i][j];
		}
	}
cerr << 11 << endl;

	// parameters
	loadPrmt();

	// adjList
	vector<int> adjList;
	int flag = 0;
	for (int i = 0; i < s.size(); i++) {
		ss += s[i];
		if (i > 0 && s[i - 1] != '&' && s[i] != '&' && !(s[i - 1] == '(' && s[i] == ')')) {
			adjList.push_back(i - 1);
		}
		if (s[i] == '(') {
			flag++;
		}
	}
cerr << 22 << endl;
	// basepairs
	bp = new int[len];
	for (int i = 0; i < len; i++) {
		bp[i] = -1;
	}
	basepairs = new Matr_(flag, 2);
	vector<int> list1;
	vector<int> list2;
	flag = 0;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '(') {
			list1.push_back(i);
		} else if (ss[i] == '[') {
			list2.push_back(i);
		} else if (ss[i] == ')') {
			basepairs->data[flag][0] = list1.back();
			basepairs->data[flag][1] = i;
			bp[i] = list1.back();
			bp[list1.back()] = i;
			flag++;
			list1.pop_back();
		} else if (ss[i] == ']') {
			basepairs->data[flag][0] = list2.back();
			basepairs->data[flag][1] = i;
			bp[i] = list2.back();
			bp[list2.back()] = i;
			flag++;
			list2.pop_back();
		}
	}
cerr << 33 << endl;
	if (t == NULL) {
		tertiary_contacts = NULL;
	} else {
		tertiary_contacts = new Matr_(t);
	}
	cycles = 1000;

	// dList
	dList.push_back(5.78);
	dList.push_back(13.7);
	if (t != NULL) {
		for (int i = 0; i < t->row; i++) {
			dList.push_back(t->data[i][2]);
		}
	}
cerr << 44 << endl;
	// bounds
	bounds = new Matr_(len, len);
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (i == j) {
				bounds->data[i][j] = 0;
			}
			bounds->data[i][j] = -1;
			bounds->data[j][i] = -1;
		}
	}
	for (int i = 0; i < adjList.size(); i++) {
		int a = adjList[i];
		int b = a + 1;
		bounds->data[a][b] = 0;
		bounds->data[b][a] = 0;
	}
	for (int i = 0; i < basepairs->row; i++) {
		int a = int(basepairs->data[i][0]);
		int b = int(basepairs->data[i][1]);
		bounds->data[a][b] = 1;
		bounds->data[b][a] = 1;
	}
	if (t != NULL) {
		for (int i = 0; i < t->row; i++) {
			int a = int(t->data[i][0]);
			int b = int(t->data[i][1]);
			bounds->data[a][b] = i + 2;
			bounds->data[b][a] = i + 2;
		}
	}
}

NAST::~NAST() {
	delete [] coord;
	delete tertiary_contacts;
	delete bounds;
	delete basepairs;
}

Point *NAST::run() {
	bounds->print();
	Point *c = new Point[len];
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < 3; j++) {
			c[i][j] = coord[i][j];
		}
	}
	double Et = energy(c);
	double Eo = Et;
	double En = Et;
	MRand mr;
	for (int i = 0; i < cycles; i++) {
		if (i % 100 == 0) {
			cout << i << '\t' << Eo << endl;
		}
		int n = int(mr.run() * len);
		double x_ = (mr.run() - 0.5) * 2;
		double y_ = (mr.run() - 0.5) * 2;
		double z_ = (mr.run() - 0.5) * 2;
		c[n][0] += x_;
		c[n][1] += y_;
		c[n][2] += z_;
		En = energy(c);
		double deltaE = En - Eo;
		double dTemp = mr.run();
		if (dTemp < exp(-50 * deltaE / Eo)) {
			Eo = En;
		} else {
			c[n][0] -= x_;
			c[n][1] -= y_;
			c[n][2] -= z_;
		}
	}
	return c;
}

void NAST::loadPrmt() {
	string lib = getenv("RNA");
	string fileName = lib + "NAST.par";
	ifstream ifile(fileName.c_str());
	string line;
	int n = 0;
	while (ifile) {
		getline(ifile, line);
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;
		vector<string> splitLine;
		tokenize(line, splitLine);
		if (splitLine.size() == 2) {
			n++;
			if (n == 1) {
				kb_all = atof(splitLine[1].c_str());
			} else if (n == 2) {
				rb_all = atof(splitLine[1].c_str());
			} else if (n == 3) {
				kb_hel = atof(splitLine[1].c_str());
			} else if (n == 4) {
				rb_hel = atof(splitLine[1].c_str());
			} else if (n == 5) {
				ka_nhel = atof(splitLine[1].c_str());
			} else if (n == 6) {
				ra_nhel = atof(splitLine[1].c_str());
			} else if (n == 7) {
				ka_hel1 = atof(splitLine[1].c_str());
			} else if (n == 8) {
				ra_hel1 = atof(splitLine[1].c_str());
			} else if (n == 9) {
				ka_hel2 = atof(splitLine[1].c_str());
			} else if (n == 10) {
				ra_hel2 = atof(splitLine[1].c_str());
			} else if (n == 11) {
				k1d_nhel = atof(splitLine[1].c_str());
			} else if (n == 12) {
				k2d_nhel = atof(splitLine[1].c_str());
			} else if (n == 13) {
				k3d_nhel = atof(splitLine[1].c_str());
			} else if (n == 14) {
				rd_nhel = atof(splitLine[1].c_str());
			} else if (n == 15) {
				k1d_hel1 = atof(splitLine[1].c_str());
			} else if (n == 16) {
				k2d_hel1 = atof(splitLine[1].c_str());
			} else if (n == 17) {
				k3d_hel1 = atof(splitLine[1].c_str());
			} else if (n == 18) {
				rd_hel1 = atof(splitLine[1].c_str());
			} else if (n == 19) {
				k1d_hel2 = atof(splitLine[1].c_str());
			} else if (n == 20) {
				k2d_hel2 = atof(splitLine[1].c_str());
			} else if (n == 21) {
				k3d_hel2 = atof(splitLine[1].c_str());
			} else if (n == 22) {
				rd_hel2 = atof(splitLine[1].c_str());
			} else if (n == 23) {
				k1d_hel3 = atof(splitLine[1].c_str());
			} else if (n == 24) {
				k2d_hel3 = atof(splitLine[1].c_str());
			} else if (n == 25) {
				k3d_hel3 = atof(splitLine[1].c_str());
			} else if (n == 26) {
				rd_hel3 = atof(splitLine[1].c_str());
			} else if (n == 27) {
				epsilon = atof(splitLine[1].c_str());
			} else if (n == 28) {
				sigma1 = atof(splitLine[1].c_str());
			} else if (n == 29) {
				sigma2 = atof(splitLine[1].c_str());
			} else if (n == 30) {
				kpt = atof(splitLine[1].c_str());
			}
		}
	}
	ifile.close();
}

double NAST::energy(Point *c) {
	double E = 0;
	for (int i = 0; i < len; i++) {
		for (int j = i + 1; j < len; j++) {
			if (bounds->data[i][j] == -1) {
				double d = c[i].dist(c[j]);
				E += 4 * epsilon * pow(sigma1 / d, 12.);
			} else if (bounds->data[i][j] == 0) {
				double d = c[i].dist(c[j]);
				E += kb_all * (d - rb_all) * (d - rb_all);
			} else if (bounds->data[i][j] == 1) {
				double d = c[i].dist(c[j]);
				E += kb_hel * (d - rb_hel) * (d - rb_hel);
			} else {
				int n = int(bounds->data[i][j]);
				double d = c[i].dist(c[j]);
				E += kpt * (d - dList[n]) * (d - dList[n]);
			}
			if (j == i + 1 && j + 1 < len) {
				if (bounds->data[i][j] == 0 && bounds->data[j][j + 1] == 0) {
					double d = Point::angle(c[i], c[j], c[j + 1]) / 180. * 3.1415927;
					if (bp[i] + i == bp[j] + j && bp[j] + j == bp[j + 1] + j + 1) {
						E += ka_hel1 * (d - ra_hel1) * (d - ra_hel1);
					} else {
						E += ka_nhel * (d - ra_nhel) * (d - ra_nhel);
					}
				}
			}
			if (j == i + 1 && bounds->data[i][j] == 0 && bp[i] != -1 && bp[j] != -1 && bp[i] + i == bp[j] + j) {
				double d = Point::angle(c[bp[i]], c[i], c[j]) / 180. * 3.1415927;
				E += ka_hel2 * (d - ra_hel2) * (d - ra_hel2);
			}

			if (j == i + 1 && j + 2 < len) {
				if (bounds->data[i][j] == 0 && bounds->data[j][j + 1] == 0 && bounds->data[j + 1][j + 2] == 0) {
					double d = (Point::dihedral(c[i], c[j], c[j + 1], c[j + 2]) - 180) / 180. * 3.1415927;
					if (bp[i] + i == bp[j] + j && bp[j] + j == bp[j + 1] + j + 1 && bp[j + 1] + j + 1 == bp[j + 2] + j + 2) {
						E += k1d_hel1 * cos(d - rd_hel1) + k2d_hel1 * cos(2 * (d - rd_hel1)) + k3d_hel1 * cos(3 * (d - rd_hel1));
					} else {
						E += k1d_nhel * cos(d - rd_nhel) + k2d_nhel * cos(2 * (d - rd_nhel)) + k3d_nhel * cos(3 * (d - rd_nhel));
					}
				}
			}
			if (j == i + 1 && bounds->data[i][j] == 0) {
				if (bp[i] != -1 && bp[j] != -1 && bp[i] + i == bp[j] + j) {
					double d = (Point::dihedral(c[bp[j]], c[bp[i]], c[i], c[j]) - 180) / 180. * 3.1415927;
					E += k1d_hel2 * cos(d - rd_hel2) + k2d_hel2 * cos(2 * (d - rd_hel2)) + k3d_hel2 * cos(3 * (d - rd_hel2));
				}
				if (j + 1 < len && bounds->data[j][j + 1] == 0 && bp[i] != -1 && bp[j] != -1 && bp[j + 1] != -1) {
					if (bp[i] + i == bp[j] + j && bp[j] + j == bp[j + 1] + j + 1 && bp[j + 1] + j + 1 == bp[j + 2] + j + 2) {
						double d = (Point::dihedral(c[bp[i]], c[i], c[j], c[j + 1]) - 180) / 180. * 3.1415927;
						E += k1d_hel3 * cos(d - rd_hel3) + k2d_hel3 * cos(2 * (d - rd_hel3)) + k3d_hel3 * cos(3 * (d - rd_hel3));
					}
				}
			}
		}
	}
	return E;
}


