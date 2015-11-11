#include "TestAtomDist.h"
#include "RNA.h"

TestAtomDist::TestAtomDist(char *str1, char *str2, int n, char *str3) {
	string string1(str1);
	string string2(str2);
	string string3(str3);
	init(string1, string2, n, string3);
}

TestAtomDist::TestAtomDist(string str1, string str2, int n, string str3) {
	init(str1, str2, n, str3);
}

void TestAtomDist::init(string name1, string name2, int n, string filename) {
	this->name1 = name1;
	this->name2 = name2;
	this->num = n;
	this->filename = filename;
	binWidth = 0.1;
	maxDist = 0;
	minDist = 9999;
	maxDist2 = 0;
	minDist2 = 9999;
	bins = int(1000 / binWidth);
	distribution = new int[bins];
	distribution2 = new int[bins];
	for (int i = 0; i < bins; i++) {
		distribution[i] = 0;
		distribution2[i] = 0;
	}
}

void TestAtomDist::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		RNA rna(str);
		int i = 0, j = 0, k = 0;
		int nj = 0;
		while (i < rna.chains.size()) {
			Atom &atom1 = rna.chains[i].residues[j].atoms[k];
			string temp1 = atom1.name;
			if (temp1 == name1) {
				int flag = num;
				int flag2 = 0;
				int i_ = i, j_ = j, nj_ = nj, k_ = k + 1;
				Residue *r1, *r2;
				while (flag > 0) {
					if (rna.chains[i_].residues[j_].atoms.size() == k_) {
						r1 = &(rna.chains[i_].residues[j_]);
						k_ = 0;
						j_++;
						nj_++;
					}
					if (rna.chains[i_].residues.size() == j_) {
						j_ = 0;
						i_++;
						if (i_ == rna.chains.size()) break;
						if (nj_ - nj == num && name2 == "P" && rna.chains[i_].residues[j_].atoms[0].name == "O5*") break;
					}
					if (k_ == 0) {
						r2 = &(rna.chains[i_].residues[j_]);
						if (!(r2->nextTo(*r1))) {
							flag2 = 1;
						}
					}

					Atom &atom2 = rna.chains[i_].residues[j_].atoms[k_];
					string temp2 = atom2.name;
					if (temp2 == name2) {
						flag--;
						if (flag == 0) {
							int n = int(atom1.dist(atom2) / binWidth);
							if (n >= 10000) break;
							if (flag2 == 0) {
								distribution[n]++;
								if (n < minDist) {
									minDist = n;
								}
								if (n > maxDist) {
									maxDist = n;
								}
							} else {
								distribution2[n]++;
								if (n < minDist2) {
									minDist2 = n;
								}
								if (n > maxDist2) {
									maxDist2 = n;
								}
							}
							break;
						}
					}

					k_++;
				}
			}
			k++;
			if (rna.chains[i].residues[j].atoms.size() == k) {
				k = 0;
				j++;
				nj++;
			}
			if (rna.chains[i].residues.size() == j) {
				j = 0;
				i++;
			}
		}
	}
	ifile.close();
	double min = minDist, max = maxDist;
	double sum = 0;
	for (int i = minDist; i <= maxDist; i++) {
		sum += distribution[i];
		cout << i * 0.1 << '\t' << distribution[i] << endl;
	}
	double sum1 = 0, sum2 = 0;
	for (int i = minDist; i <= maxDist; i++) {
		sum1 += distribution[i] / sum;
		if (sum1 > 0.025) {
			min = i * 0.1;
			break;
		}
	}
	for (int i = maxDist; i >= minDist; i--) {
		sum2 += distribution[i] / sum;
		if (sum2 > 0.025) {
			max = i * 0.1;
			break;
		}
	}
	cout << endl;

	double min2 = minDist2, max2 = maxDist2;
	sum = 0;
	for (int i = minDist2; i <= maxDist2; i++) {
		sum += distribution2[i];
		cout << i * 0.1 << '\t' << distribution2[i] << endl;
	}
	sum1 = 0; sum2 = 0;
	for (int i = minDist2; i <= maxDist2; i++) {
		sum1 += distribution2[i] / sum;
		if (sum1 > 0.025) {
			min2 = i * 0.1;
			break;
		}
	}
	for (int i = maxDist2; i >= minDist2; i--) {
		sum2 += distribution2[i] / sum;
		if (sum2 > 0.025) {
			max2 = i * 0.1;
			break;
		}
	}
	cout << endl;
	cerr << name1 << '-' << num << '-' << name2 << ' ' << min << ' ' << max << ' ' << min2 << ' ' << max2 << endl;

}


