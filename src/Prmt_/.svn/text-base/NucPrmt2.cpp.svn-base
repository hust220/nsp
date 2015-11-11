#include "NucPrmt2.h"
#include "RNA.h"

NucPrmt2::NucPrmt2(char *str3) {
	string string3(str3);
	init(string3);
}

NucPrmt2::NucPrmt2(string str3) {
	init(str3);
}

void NucPrmt2::init(string filename) {
	this->filename = filename;
	binWidth = 0.1;
	bins = int(1000 / binWidth);

	prmt = new Matr_ *[4];
	prmtMinMax = new Matr_ *[4];
	for (int i = 0; i < 4; i++) {
		prmt[i] = new Matr_(10, bins);
		prmtMinMax[i] = new Matr_(10, 2);
		for (int j = 0; j < 10; j++) {
			for (int k = 0; k < bins; k++) {
				prmt[i]->data[j][k] = 0;
			}
			prmtMinMax[i]->data[j][0] = 999;
			prmtMinMax[i]->data[j][1] = 0;
		}
	}
}

void NucPrmt2::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		cerr << str << endl;
		RNA rna(str);
		int len = rna.getLen();
		int *type = new int[len];
	
		Point **pn = new Point *[5];
		for (int i = 0; i < 5; i++) {
			pn[i] = NULL;
		}
		int i = 0, j = 0, k = 0, nj = 0;
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
				pn[0] = atom.coord();
			} else if (atom.name == "C4*") {
				pn[1] = atom.coord();
			}
			if (type[nj] == 0) {
				if (atom.name == "C4") {
					pn[2] = atom.coord();
				} else if (atom.name == "N1") {
					pn[3] = atom.coord();
				} else if (atom.name == "N6") {
					pn[4] = atom.coord();
				}
			} else if (type[nj] == 1) {
				if (atom.name == "C6") {
					pn[2] = atom.coord();
				} else if (atom.name == "N3") {
					pn[3] = atom.coord();
				} else if (atom.name == "O4") {
					pn[4] = atom.coord();
				}
			} else if (type[nj] == 2) {
				if (atom.name == "C4") {
					pn[2] = atom.coord();
				} else if (atom.name == "N1") {
					pn[3] = atom.coord();
				} else if (atom.name == "O6") {
					pn[4] = atom.coord();
				}
			} else if (type[nj] == 3) {
				if (atom.name == "C6") {
					pn[2] = atom.coord();
				} else if (atom.name == "N3") {
					pn[3] = atom.coord();
				} else if (atom.name == "N4") {
					pn[4] = atom.coord();
				}
			}
			if (k + 1 == rna[i][j].atoms.size()) {
				int iType = type[nj];
				for (int ii = 0; ii < 5; ii++) {
					if (pn[ii] == NULL) continue;
					for (int jj = ii + 1; jj < 5; jj++) {
						if (pn[jj] == NULL) continue;
						double dTemp = pn[ii]->dist(pn[jj]);
						if (dTemp > 1000) continue;
						int iTemp = int(dTemp / binWidth);
						int n = ii * 5 + jj - (ii + 1) * (ii + 2) / 2;
						if (iTemp < prmtMinMax[iType]->data[n][0]) {
							prmtMinMax[iType]->data[n][0] = iTemp;
						}
						if (iTemp > prmtMinMax[iType]->data[n][1]) {
							prmtMinMax[iType]->data[n][1] = iTemp;
						}
						prmt[iType]->data[n][iTemp]++;
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

				for (int i = 0; i < 5; i++) {
					delete pn[i];
				}
				delete [] pn;
				pn = new Point *[5];
				for (int i = 0; i < 5; i++) {
					pn[i] = NULL;
				}
			}
		}
	}
	ifile.close();

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 10; j++) {
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

