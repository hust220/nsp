#include "BPPrmt.h"
#include "RNA.h"

BPPrmt::BPPrmt(char *str1) {
	string string1(str1);
	init(string1);
}

BPPrmt::BPPrmt(string str1) {
	init(str1);
}

void BPPrmt::init(string filename) {
	this->filename = filename;
	binWidth = 0.1;
	bins = int(1000 / binWidth);
	prmt = new Matr_ *[16];
	prmtMinMax = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		prmt[i] = new Matr_(169, bins);
		prmtMinMax[i] = new Matr_(169, 2);
		for (int j = 0; j < 169; j++) {
			for (int k = 0; k < bins; k++) {
				prmt[i]->data[j][k] = 0;
			}
			prmtMinMax[i]->data[j][0] = 999;
			prmtMinMax[i]->data[j][1] = 0;
		}
	}
}

void BPPrmt::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		RNA rna(str);
		cerr << str << endl;
		int len = rna.getLen();
		if (len < 6) continue;
		Point **p1 = new Point *[13];
		Point **p2 = new Point *[13];
		for (int i = 0; i < 13; i++) {
			p1[i] = NULL;
			p2[i] = NULL;
		}
		int type1, type2;
		int i = 0, j = 0, k = 0;
		int nj = 0;
		while (i < rna.chains.size()) {
			Atom &atom = rna[i][j][k];
			if (nj == 1) {
				if (rna[i][j][0].name == "P" && k < 13) {
					p1[k] = atom.coord();
				} else if (rna[i][j][0].name != "P" && k < 10) {
					p1[k + 3] = atom.coord();
				}
			}
			if (nj + 2 == len) {
				if (rna[i][j][0].name == "P" && k < 13) {
					p2[k] = atom.coord();
				} else if (rna[i][j][0].name != "P" && k < 10) {
					p2[k + 3] = atom.coord();
				}
			}

			k++;
			if (k == rna.chains[i].residues[j].atoms.size()) {
				k = 0;
				j++;
				nj++;
				if (j == rna.chains[i].residues.size()) {
					j = 0;
					i++;
					if (i == rna.chains.size()) break;
				}
				if (nj == 1) {
					if (rna[i][j].name == "A") {
						type1 = 0;
					} else if (rna[i][j].name == "U") {
						type1 = 1;
					} else if (rna[i][j].name == "G") {
						type1 = 2;
					} else {
						type1 = 3;
					}
				}
				if (nj == len - 2) {
					if (rna[i][j].name == "A") {
						type2 = 0;
					} else if (rna[i][j].name == "U") {
						type2 = 1;
					} else if (rna[i][j].name == "G") {
						type2 = 2;
					} else {
						type2 = 3;
					}
				}
			}
		}

		int type = type1 * 4 + type2;
		for (int ii = 0; ii < 13; ii++) {
			if (p1[ii] == NULL) continue;
			for (int jj = 0; jj < 13; jj++) {
				if (p2[jj] == NULL) continue;
				double dTemp = p1[ii]->dist(p2[jj]);
				if (dTemp > 999) continue;
				int iTemp = int(dTemp / binWidth);
				if (iTemp < prmtMinMax[type]->data[ii * 13 + jj][0]) {
					prmtMinMax[type]->data[ii * 13 + jj][0] = iTemp;
				}
				if (iTemp > prmtMinMax[type]->data[ii * 13 + jj][1]) {
					prmtMinMax[type]->data[ii * 13 + jj][1] = iTemp;
				}
				prmt[type]->data[ii * 13 + jj][iTemp]++;
			}
		}
		for (int i = 0; i < 13; i++) {
			delete p1[i];
			delete p2[i];
		}
		delete [] p1;
		delete [] p2;
	}
	ifile.close();

	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 13 * 13; j++) {
			double min = prmtMinMax[i]->data[j][0];
			double max = prmtMinMax[i]->data[j][1];
			double sum = 0;
			double sum1 = 0, sum2 = 0;
			for (int k = int(prmtMinMax[i]->data[j][0]); k <= int(prmtMinMax[i]->data[j][1]); k++) {
				sum += prmt[i]->data[j][k];
			}
			for (int k = int(prmtMinMax[i]->data[j][0]); k <= int(prmtMinMax[i]->data[j][1]); k++) {
				sum1 += prmt[i]->data[j][k] / sum;
				if (sum1 > 0.025) {
					min = k * binWidth;
					break;
				}
			}
			for (int k = int(prmtMinMax[i]->data[j][1]); k >= int(prmtMinMax[i]->data[j][0]); k--) {
				sum2 += prmt[i]->data[j][k] / sum;
				if (sum2 > 0.025) {
					max = k * binWidth;
					break;
				}
			}
		//	prmtMinMax[i]->data[j][0] = min;
		//	prmtMinMax[i]->data[j][1] = max;
			cout << min << ' ' << max << endl;
		}
	}
}


