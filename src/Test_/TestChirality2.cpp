#include "TestChirality2.h"
#include "RNA.h"

TestChirality2::TestChirality2(char *str3) {
	string string3(str3);
	init(string3);
}

TestChirality2::TestChirality2(string str3) {
	init(str3);
}

void TestChirality2::init(string filename) {
	this->filename = filename;
	minDist = 199;
	maxDist = 0;
	binWidth = 0.1;
	bins = int(200 / binWidth);
	chir = new double[bins];
	for (int i = 0; i < bins; i++) {
		chir[i] = 0;
	}
}

void TestChirality2::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		cerr << str << endl;
		Point *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL;
		string s[4];
		RNA rna(str);
		int i = 0, j = 0, k = 0;
		int nj = 0;
		while (i < rna.chains.size()) {
			Atom &atom = rna.chains[i].residues[j].atoms[k];
			if (atom.name == "C3*") {
				delete p1;
				p1 = atom.coord();
			} else if (atom.name == "C2*") {
				delete p2;
				p2 = atom.coord();
			} else if (atom.name == "O2*") {
				p3 = atom.coord();
			} else if (atom.name == "C1*") {
				p4 = atom.coord();
			}
			if (k + 1 == rna[i][j].atoms.size()) {
				int c = int((Point::chirality(p1, p2, p3, p4) + 100) / binWidth);
				chir[c]++;
				if (c < minDist) {
					minDist = c;
				}
				if (c > maxDist) {
					maxDist = c;
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
		delete p1;
		p1 = NULL;
		delete p2;
		p2 = NULL;
		delete p3;
		p3 = NULL;
		delete p4;
		p4 = NULL;
	}
	ifile.close();

	double min = minDist;
	double max = maxDist;
	double center;
	double sum = 0;
	double temp = 0;
	for (int i = minDist; i <= maxDist; i++) {
		if (temp < chir[i]) {
			temp = chir[i];
			center = (i - 1000) * 0.1;
		}
		sum += chir[i];
//		cout << (i - 1000) * 0.1 << '\t' << chir->data[n][i] << endl;
	}
	double sum1 = 0, sum2 = 0;
	for (int i = minDist; i <= maxDist; i++) {
		sum1 += chir[i] / sum;
		if (sum1 > 0.025) {
			min = (i - 1000) * 0.1;
			break;
		}
	}
	for (int i = maxDist; i >= minDist; i--) {
		sum2 += chir[i] / sum;
		if (sum2 > 0.025) {
			max = (i - 1000) * 0.1;
			break;
		}
	}
	for (int i = int(min * 10 + 1000); i <= int(max * 10 + 1000); i++) {
		cout << (i - 1000) * 0.1 << ' ' << chir[i] << endl;
	}
	cout << endl;

	cerr << min << ' ' << max << ' ' << center << endl;
}


