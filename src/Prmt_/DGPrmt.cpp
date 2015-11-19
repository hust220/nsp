#include "DGPrmt.h"
#include "RNA.h"

DGPrmt::DGPrmt(char *str3) {
	string string3(str3);
	init(string3);
}

DGPrmt::DGPrmt(string str3) {
	init(str3);
}

void DGPrmt::init(string filename) {
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

void DGPrmt::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		RNA rna(str);
		cerr << str << endl;
		int len = rna.getLen();
		int type[len];
		for (int i = 0, flag = 0; i < rna.chains.size(); i++) {
			for (int j = 0; j < rna[i].residues.size(); j++, flag++) {
				if (rna[i][j].name == "A") {
					type[flag] = 0;
				} else if (rna[i][j].name == "U") {
					type[flag] = 1;
				} else if (rna[i][j].name == "G") {
					type[flag] = 2;
				} else {
					type[flag] = 3;
				}
			}
		}
		int i = 0, j = 0, nj = 0, k = 0;
		Residue *rn = &(rna[0][0]), *ro = NULL;
		Point **pn = new Point *[13];
		Point **po = new Point *[13];
		for (int i = 0; i < 13; i++) {
			pn[i] = NULL;
			po[i] = NULL;
		}
		while (i < rna.chains.size()) {
			Atom &atom = rna[i][j][k];
			if (rna[i][j][0].name != "P" && k < 10) {
				pn[3 + k] = atom.coord();
			} else if (k < 13) {
				pn[k] = atom.coord();
			}
			if (k + 1 == rna[i][j].atoms.size() && nj != 0) {
				if (rn->nextTo(*ro)) {
					int iType = type[nj - 1] * 4 + type[nj];
					for (int ii = 0; ii < 13; ii++) {
						if (po[ii] == NULL) continue;
						for (int jj = 0; jj < 13; jj++) {
							if (pn[jj] == NULL) continue;
							double dTemp = po[ii]->dist(pn[jj]);
							if (dTemp > 1000) continue;
							int iTemp = int(dTemp / binWidth);
							if (iTemp < prmtMinMax[iType]->data[ii * 13 + jj][0]) {
								prmtMinMax[iType]->data[ii * 13 + jj][0] = iTemp;
							}
							if (iTemp > prmtMinMax[iType]->data[ii * 13 + jj][1]) {
								prmtMinMax[iType]->data[ii * 13 + jj][1] = iTemp;
							}
							prmt[iType]->data[ii * 13 + jj][iTemp]++;
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

				ro = rn;
				rn = &(rna[i][j]);
				for (int i = 0; i < 13; i++) {
					delete po[i];
				}
				delete [] po;
				po = pn;
				pn = new Point *[13];
				for (int i = 0; i < 13; i++) {
					pn[i] = NULL;
				}
			}
		}
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


