#include "TestChirality.h"
#include "RNA.h"

TestChirality::TestChirality(char *str3) {
	string string3(str3);
	init(string3);
}

TestChirality::TestChirality(string str3) {
	init(str3);
}

void TestChirality::init(string filename) {
	this->filename = filename;
	for (int i = 0; i < 6; i++) {
		minDist[i] = 199;
		maxDist[i] = 0;
	}
	binWidth = 0.1;
	bins = int(200 / binWidth);
	chir = new Matr_(6, bins);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < bins; j++) {
			chir->data[i][j] = 0;
		}
	}
}

void TestChirality::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		Point *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL;
		string s[4];
		RNA rna(str);
		int i = 0, j = 0, k = 0;
		int nj = 0;
		while (i < rna.chains.size()) {
			Atom &atom = rna.chains[i].residues[j].atoms[k];
			if (atom.name == "P" || atom.name == "O5*" || atom.name == "C5*" || atom.name == "C4*" || atom.name == "C3*" || atom.name == "O3*") {
				delete p1;
				p1 = p2;
				s[0] = s[1];
				p2 = p3;
				s[1] = s[2];
				p3 = p4;
				s[2] = s[3];
				p4 = atom.coord();
				s[3] = atom.name;
			}
			if (atom.name == "P" && p3 != NULL) {
				if (p4->dist(p3) > 2) {
					delete p1;
					p1 = NULL;
					s[0] = "";
					delete p2;
					p2 = NULL;
					s[1] = "";
					delete p3;
					p3 = NULL;
					s[2] = "";
				}
			}
			if (atom.name == "O5*") {
				if (s[2] == "O3*") {
					delete p1;
					p1 = NULL;
					s[0] = "";
					delete p2;
					p2 = NULL;
					s[1] = "";
					delete p3;
					p3 = NULL;
					s[2] = "";
				}
			}
			if (p1 != NULL && p2 != NULL && p3 != NULL && p4 != NULL) {
				int flag;
				if (s[1] == "P") {
					flag = 0;
				} else if (s[1] == "O5*") {
					flag = 1;
				} else if (s[1] == "C5*") {
					flag = 2;
				} else if (s[1] == "C4*") {
					flag = 3;
				} else if (s[1] == "C3*") {
					flag = 4;
				} else if (s[1] == "O3*") {
					flag = 5;
				}
				int c = int((Point::chirality(p1, p2, p3, p4) + 100) / binWidth);
				chir->data[flag][c]++;
				if (c < minDist[flag]) {
					minDist[flag] = c;
				}
				if (c > maxDist[flag]) {
					maxDist[flag] = c;
				}
			}
			k++;
			if (rna.chains[i].residues[j].atoms.size() == k) {
				k = 0;
				j++;
			}
			if (rna.chains[i].residues.size() == j) {
				j = 0;
				i++;
			}
		}
		delete p1;
		delete p2;
		delete p3;
		delete p4;
	}
	ifile.close();

	double min[6];
	double max[6];
	double center[6];
	for (int n = 0; n < 6; n++) {
		min[n] = minDist[n];
		max[n] = maxDist[n];
		double sum = 0;
		double temp = 0;
		for (int i = minDist[n]; i <= maxDist[n]; i++) {
			if (temp < chir->data[n][i]) {
				temp = chir->data[n][i];
				center[n] = (i - 1000) * 0.1;
			}
			sum += chir->data[n][i];
//			cout << (i - 1000) * 0.1 << '\t' << chir->data[n][i] << endl;
		}
		double sum1 = 0, sum2 = 0;
		for (int i = minDist[n]; i <= maxDist[n]; i++) {
			sum1 += chir->data[n][i] / sum;
			if (sum1 > 0.025) {
				min[n] = (i - 1000) * 0.1;
				break;
			}
		}
		for (int i = maxDist[n]; i >= minDist[n]; i--) {
			sum2 += chir->data[n][i] / sum;
			if (sum2 > 0.025) {
				max[n] = (i - 1000) * 0.1;
				break;
			}
		}
		for (int i = int(min[n] * 10 + 1000); i <= int(max[n] * 10 + 1000); i++) {
			cout << (i - 1000) * 0.1 << ' ' << chir->data[n][i] << endl;
		}
		cout << endl;
	}
	for (int i = 0; i < 6; i++) {
		cerr << min[i] << ' ' << max[i] << ' ' << center[i] << endl;
	}

}


