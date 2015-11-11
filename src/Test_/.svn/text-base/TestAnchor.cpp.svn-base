#include "TestAnchor.h"
#include "RNA.h"

TestAnchor::TestAnchor(char *str1, char *str2, char *str3) {
	string string1(str1);
	string string2(str2);
	string string3(str3);
	init(string1, string2, string3);
}

TestAnchor::TestAnchor(string str1, string str2, string str3) {
	init(str1, str2, str3);
}

void TestAnchor::init(string name1, string name2, string filename) {
	this->name1 = name1;
	this->name2 = name2;
	this->filename = filename;
	binWidth = 0.1;
	maxDist = 0;
	minDist = 9999;
	bins = int(1000 / binWidth);
	distribution = new int[bins];
	for (int i = 0; i < bins; i++) {
		distribution[i] = 0;
	}
}

void TestAnchor::run() {
	ifstream ifile(filename.c_str());
	string str;
	while (ifile >> str) {
		RNA rna(str);
		int len = rna.getLen();
		int i = 0, j = 0, k = 0;
		int nj = 0;
		Point *p1 = NULL, *p2 = NULL;
		while (i < rna.chains.size()) {
			Atom &atom = rna.chains[i].residues[j].atoms[k];
			if (nj == 0 && atom.name == name1) {
				p1 = atom.coord();
			}
			if (nj == len - 1 && atom.name == name2) {
				p2 = atom.coord();
			}
			k++;
			if (k == rna.chains[i].residues[j].atoms.size()) {
				k = 0;
				j++;
				nj++;
			}
			if (j == rna.chains[i].residues.size()) {
				j = 0;
				i++;
			}
		}
		if (p1 == NULL || p2 == NULL) {
			delete p1;
			delete p2;
			continue;
		}
		double distance = p1->dist(p2);
		delete p1;
		delete p2;
		int n = int(distance / binWidth);
		if (n >= 10000) continue;
		distribution[n]++;
		if (n < minDist) {
			minDist = n;
		}
		if (n > maxDist) {
			maxDist = n;
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

	cerr << name1 << '-' << name2 << ' ' << min << ' ' << max << endl;

}


