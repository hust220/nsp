#include "NucPrmt3.h"
#include "RNA.h"

NucPrmt3::NucPrmt3(char *str3) {
	string string3(str3);
	init(string3);
}

NucPrmt3::NucPrmt3(string str3) {
	init(str3);
}

void NucPrmt3::init(string filename) {
	this->filename = filename;
	binWidth = 0.1;
	bins = int(1000 / binWidth);

	prmt = new Matr_ *[4];
	prmtMinMax = new Matr_ *[4];
	lens = new int[4];
	lens[0] = 23;
	lens[1] = 21;
	lens[2] = 24;
	lens[3] = 21;

	for (int i = 0; i < 4; i++) {
		int l = lens[i] * (lens[i] - 1) / 2;
		prmt[i] = new Matr_(l, bins);
		prmtMinMax[i] = new Matr_(l, 2);
		for (int j = 0; j < l; j++) {
			for (int k = 0; k < bins; k++) {
				prmt[i]->data[j][k] = 0;
			}
			prmtMinMax[i]->data[j][0] = 999;
			prmtMinMax[i]->data[j][1] = 0;
		}
	}
}

void NucPrmt3::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		cerr << str << endl;
		RNA rna(str);
		int len = rna.getLen();
		int *type = new int[len];
	
		int i = 0, j = 0, k = 0, nj = 0;
		Point **pn = NULL;
		while (i < rna.chains.size()) {
			Atom &atom = rna[i][j][k];
			if (k == 0) {
				if (nj > 0) {
					int length = lens[type[nj - 1]];
					pn[length - 1] = atom.coord();
					int iType = type[nj - 1];
					for (int ii = 0; ii < length; ii++) {
						if (pn[ii] == NULL) continue;
						for (int jj = ii + 1; jj < length; jj++) {
							if (pn[jj] == NULL) continue;
							double dTemp = pn[ii]->dist(pn[jj]);
							if (dTemp > 1000) continue;
							int iTemp = int(dTemp / binWidth);
							int n = ii * length + jj - (ii + 1) * (ii + 2) / 2;
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
				if (rna[i][j].name == "A") {
					type[nj] = 0;
				} else if (rna[i][j].name == "U") {
					type[nj] = 1;
				} else if (rna[i][j].name == "G") {
					type[nj] = 2;
				} else {
					type[nj] = 3;
				}
				if (pn != NULL) {
					for (int i = 0; i < lens[type[nj - 1]]; i++) {
						delete pn[i];
					}
					delete [] pn;
					if (nj == len - 1) break;
				}
				pn = new Point *[lens[type[nj]]];
				for (int i = 0; i < lens[type[nj]]; i++) {
					pn[i] = NULL;
				}
			}
			if (rna[i][j][0].name == "P") {
				pn[k] = atom.coord();
			} else if (rna[i][j][0].name == "C4*") {
				pn[k + 3] = atom.coord();
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

			}
		}
	}
	ifile.close();

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < lens[i] * (lens[i] - 1) / 2.0; j++) {
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

