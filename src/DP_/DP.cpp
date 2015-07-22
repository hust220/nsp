#include "DP.h"

DP::DP(string seq) {
	max_span_ = 100;
	min_hairpin_len_ = 5;

	seq_ = seq;
	len_ = seq.size();

	types_ = new int[len_];

	int temp1 = len_ - max_span_ + 1;
	int temp2 = len_ - min_hairpin_len_ + 1;
	int temp = (temp1 + temp2) * (temp2 - temp1 + 1) / 2;
	en_ = new double[temp];

	f_ = new double[len_];
}

void DP::run() {
}

void DP::setF() {
}

void DP::setEn(int m, int n) {
	setTypes();

	int temp1 = max_span_ - min_hairpin_len_ + 1;
	int temp2 = min_hairpin_len_ - 1;

	for (int i = 0; i < temp1; i++) {
		for (int j = 0; j < len_ - min_hairpin_len_ + 1; j++) {
			int m = j;
			int n = j + i + temp2;
			double min;
			for (int ii = 0; ii < n - m + 1; ii++) {
				double temp = getEn(m, m + ii - 1) + getEn(m + ii + 1, n - 1) + getPairEn(types_[m + ii], types_[n]);
				if (min > temp) {
					min = temp;
				}
			}
			en_[m * temp1 + (n - (j + temp2))] = min;
		}
	}
}

double DP::getEn(int m, int n) {
	if (n - m + 1 < min_hairpin_len_) {
		return 0;
	}
	if (n - m + 1 > max_span_) {
		cerr << "DP::getEn error!" << endl;
		exit(1);
	}

	int temp1 = max_span_ - min_hairpin_len_ + 1;
	int temp2 = min_hairpin_len_ - 1;
	return en_[m * temp1 + (n - (m + temp2))];
}



