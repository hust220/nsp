#include "NxtPrmt2.h"
#include "RNA.h"

NxtPrmt2::NxtPrmt2(char *str3) {
	string string3(str3);
	init(string3);
}

NxtPrmt2::NxtPrmt2(string str3) {
	init(str3);
}

void NxtPrmt2::init(string filename) {
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

void NxtPrmt2::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		RNA rna(str);
		cerr << str << endl;
		int len = rna.getLen();
		int type[len];
		int i = 0, j = 0, nj = 0, k = 0;
		Residue *r1 = NULL;
		Residue *r2 = NULL;
		Residue *r3 = &(rna[0][0]);
		Point **p1 = new Point *[5];
		Point **p2 = new Point *[5];
		Point **p3 = new Point *[5];
		for (int i = 0; i < 5; i++) {
			p1[i] = NULL;
			p2[i] = NULL;
			p3[i] = NULL;
		}
		while (i < rna.chains.size()) {
			Atom &atom = rna[i][j][k];
			if (k == 0) {
				if (rna[i][j].name == "A") {
					type[nj] = 0;
				} else if (rna[i][j].name == "U") {
					type[nj] = 1;
				} else if (rna[i][j].name == "G") {
					type[nj] = 2;
				} else {
					type[nj] = 3;
				}
			}
			if (atom.name == "P") {
				p3[0] = atom.coord();
			} else if (atom.name == "C4*") {
				p3[1] = atom.coord();
			}
			if (type[nj] == 0) {
				if (atom.name == "C4") {
					p3[2] = atom.coord();
				} else if (atom.name == "N1") {
					p3[3] = atom.coord();
				} else if (atom.name == "N6") {
					p3[4] = atom.coord();
				}
			} else if (type[nj] == 1) {
				if (atom.name == "C6") {
					p3[2] = atom.coord();
				} else if (atom.name == "N3") {
					p3[3] = atom.coord();
				} else if (atom.name == "O4") {
					p3[4] = atom.coord();
				}
			} else if (type[nj] == 2) {
				if (atom.name == "C4") {
					p3[2] = atom.coord();
				} else if (atom.name == "N1") {
					p3[3] = atom.coord();
				} else if (atom.name == "O6") {
					p3[4] = atom.coord();
				}
			} else if (type[nj] == 3) {
				if (atom.name == "C6") {
					p3[2] = atom.coord();
				} else if (atom.name == "N3") {
					p3[3] = atom.coord();
				} else if (atom.name == "N4") {
					p3[4] = atom.coord();
				}
			}

			if (k + 1 == rna[i][j].atoms.size() && nj > 1) {
				if (r3->nextTo(*r2) && r2->nextTo(*r1)) {
					int iType = type[nj - 1] * 4 + type[nj];
					for (int ii = 0; ii < 5; ii++) {
						if (p1[ii] == NULL) continue;
						for (int jj = 0; jj < 5; jj++) {
							if (p3[jj] == NULL) continue;
							double dTemp = p1[ii]->dist(p3[jj]);
							if (dTemp > 1000) continue;
							int iTemp = int(dTemp / binWidth);
							if (iTemp < prmtMinMax[iType]->data[ii * 5 + jj][0]) {
								prmtMinMax[iType]->data[ii * 5 + jj][0] = iTemp;
							}
							if (iTemp > prmtMinMax[iType]->data[ii * 5 + jj][1]) {
								prmtMinMax[iType]->data[ii * 5 + jj][1] = iTemp;
							}
							prmt[iType]->data[ii * 5 + jj][iTemp]++;
						}
					}
				}
			}

			k++;
			if (rna[i][j].atoms.size() == k) {
				k = 0;
				j++;
				nj++;
				if (rna[i].residues.size() == j) {
					j = 0;
					i++;
					if (i == rna.chains.size()) break;
				}

				r1 = r2;
				r2 = r3;
				r3 = &(rna[i][j]);
				for (int i = 0; i < 5; i++) {
					delete p1[i];
				}
				delete [] p1;
				p1 = p2;
				p2 = p3;
				p3 = new Point *[5];
				for (int i = 0; i < 5; i++) {
					p3[i] = NULL;
				}
			}
		}
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


