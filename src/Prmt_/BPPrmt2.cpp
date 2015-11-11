#include "BPPrmt2.h"
#include "RNA.h"

BPPrmt2::BPPrmt2(char *str1) {
	string string1(str1);
	init(string1);
}

BPPrmt2::BPPrmt2(string str1) {
	init(str1);
}

void BPPrmt2::init(string filename) {
	this->filename = filename;
	binWidth = 0.1;
	bins = int(1000 / binWidth);
	prmt = new Matr_ *[16];
	prmtMinMax = new Matr_ *[16];
	for (int i = 0; i < 16; i++) {
		prmt[i] = new Matr_(25, bins);
		prmtMinMax[i] = new Matr_(25, 2);
		for (int j = 0; j < 25; j++) {
			for (int k = 0; k < bins; k++) {
				prmt[i]->data[j][k] = 0;
			}
			prmtMinMax[i]->data[j][0] = 999;
			prmtMinMax[i]->data[j][1] = 0;
		}
	}
}

void BPPrmt2::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		RNA rna(str);
		cerr << str << endl;
		int len = rna.getLen();
		if (len < 6) continue;
		Point **p1 = new Point *[5];
		Point **p2 = new Point *[5];
		for (int i = 0; i < 5; i++) {
			p1[i] = NULL;
			p2[i] = NULL;
		}
		int type1, type2;
		int i = 0, j = 0, k = 0;
		int nj = 0;
		while (i < rna.chains.size()) {
			Atom &atom = rna[i][j][k];
			if (nj == 1 && k == 0) {
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
			if (nj == 1) {
				if (atom.name == "P") {
					p1[0] = atom.coord();
				} else if (atom.name == "C4*") {
					p1[1] = atom.coord();
				}
				if (type1 == 0) {
					if (atom.name == "C4") {
						p1[2] = atom.coord();
					} else if (atom.name == "N1") {
						p1[3] = atom.coord();
					} else if (atom.name == "N6") {
						p1[4] = atom.coord();
					}
				} else if (type1 == 1) {
					if (atom.name == "C6") {
						p1[2] = atom.coord();
					} else if (atom.name == "N3") {
						p1[3] = atom.coord();
					} else if (atom.name == "O4") {
						p1[4] = atom.coord();
					}
				} else if (type1 == 2) {
					if (atom.name == "C4") {
						p1[2] = atom.coord();
					} else if (atom.name == "N1") {
						p1[3] = atom.coord();
					} else if (atom.name == "O6") {
						p1[4] = atom.coord();
					}
				} else if (type1 == 3) {
					if (atom.name == "C6") {
						p1[2] = atom.coord();
					} else if (atom.name == "N3") {
						p1[3] = atom.coord();
					} else if (atom.name == "N4") {
						p1[4] = atom.coord();
					}
				}
			}
			if (nj == len - 2 && k == 0) {
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
			if (nj + 2 == len) {
				if (atom.name == "P") {
					p2[0] = atom.coord();
				} else if (atom.name == "C4*") {
					p2[1] = atom.coord();
				}
				if (type2 == 0) {
					if (atom.name == "C4") {
						p2[2] = atom.coord();
					} else if (atom.name == "N1") {
						p2[3] = atom.coord();
					} else if (atom.name == "N6") {
						p2[4] = atom.coord();
					}
				} else if (type2 == 1) {
					if (atom.name == "C6") {
						p2[2] = atom.coord();
					} else if (atom.name == "N3") {
						p2[3] = atom.coord();
					} else if (atom.name == "O4") {
						p2[4] = atom.coord();
					}
				} else if (type2 == 2) {
					if (atom.name == "C4") {
						p2[2] = atom.coord();
					} else if (atom.name == "N1") {
						p2[3] = atom.coord();
					} else if (atom.name == "O6") {
						p2[4] = atom.coord();
					}
				} else if (type2 == 3) {
					if (atom.name == "C6") {
						p2[2] = atom.coord();
					} else if (atom.name == "N3") {
						p2[3] = atom.coord();
					} else if (atom.name == "N4") {
						p2[4] = atom.coord();
					}
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
			}
		}

		int type = type1 * 4 + type2;
		for (int ii = 0; ii < 5; ii++) {
			if (p1[ii] == NULL) continue;
			for (int jj = 0; jj < 5; jj++) {
				if (p2[jj] == NULL) continue;
				double dTemp = p1[ii]->dist(p2[jj]);
				if (dTemp > 999) continue;
				int iTemp = int(dTemp / binWidth);
				if (iTemp < prmtMinMax[type]->data[ii * 5 + jj][0]) {
					prmtMinMax[type]->data[ii * 5 + jj][0] = iTemp;
				}
				if (iTemp > prmtMinMax[type]->data[ii * 5 + jj][1]) {
					prmtMinMax[type]->data[ii * 5 + jj][1] = iTemp;
				}
				prmt[type]->data[ii * 5 + jj][iTemp]++;
			}
		}
		for (int i = 0; i < 5; i++) {
			delete p1[i];
			delete p2[i];
		}
		delete [] p1;
		delete [] p2;
	}
	ifile.close();

	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 5 * 5; j++) {
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


