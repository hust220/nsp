#include "NucPrmt4.h"
#include "RNA.h"

NucPrmt4::NucPrmt4(char *str3) {
	string string3(str3);
	init(string3);
}

NucPrmt4::NucPrmt4(string str3) {
	init(str3);
}

void NucPrmt4::init(string filename) {
	this->filename = filename;
	binWidth = 0.1;
	bins = int(1000 / binWidth);

	int l = 14 * (14 - 1) / 2;
	prmt = new Matr_(l, bins);
	prmtMinMax = new Matr_(l, 2);
	for (int j = 0; j < l; j++) {
		for (int k = 0; k < bins; k++) {
			prmt->data[j][k] = 0;
		}
		prmtMinMax->data[j][0] = 999;
		prmtMinMax->data[j][1] = 0;
	}
}

void NucPrmt4::run() {
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
					pn[13] = atom.coord();
					for (int ii = 0; ii < 14; ii++) {
						if (pn[ii] == NULL) continue;
						for (int jj = ii + 1; jj < 14; jj++) {
							if (pn[jj] == NULL) continue;
							double dTemp = pn[ii]->dist(pn[jj]);
							if (dTemp > 1000) continue;
							int iTemp = int(dTemp / binWidth);
							int n = ii * 14 + jj - (ii + 1) * (ii + 2) / 2;
							if (iTemp < prmtMinMax->data[n][0]) {
								prmtMinMax->data[n][0] = iTemp;
							}
							if (iTemp > prmtMinMax->data[n][1]) {
								prmtMinMax->data[n][1] = iTemp;
							}
							prmt->data[n][iTemp]++;
						}
					}
				}
				if (pn != NULL) {
					for (int i = 0; i < 14; i++) {
						delete pn[i];
					}
					delete [] pn;
					if (nj == len - 1) break;
				}
				pn = new Point *[14];
				for (int i = 0; i < 14; i++) {
					pn[i] = NULL;
				}
			}
			if (rna[i][j][0].name == "P" && k < 13) {
				pn[k] = atom.coord();
			} else if (rna[i][j][0].name == "C4*" && k < 10) {
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

	for (int j = 0; j < 14 * (14 - 1) / 2.0; j++) {
		double min = prmtMinMax->data[j][0];
		double max = prmtMinMax->data[j][1];
		double sum = 0;
		double sum1 = 0, sum2 = 0;
		for (int k = int(prmtMinMax->data[j][0]); k <= int(prmtMinMax->data[j][1]); k++) {
			sum += prmt->data[j][k];
		}
		for (int k = int(prmtMinMax->data[j][0]); k <= int(prmtMinMax->data[j][1]); k++) {
			sum1 += prmt->data[j][k] / sum;
			if (sum1 > 0.025) {
				min = k * binWidth;
				break;
			}
		}
		for (int k = int(prmtMinMax->data[j][1]); k >= int(prmtMinMax->data[j][0]); k--) {
			sum2 += prmt->data[j][k] / sum;
			if (sum2 > 0.025) {
				max = k * binWidth;
				break;
			}
		}
		cout << min << ' ' << max << endl;
	}
}

