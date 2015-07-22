#include "DistAnal.h"

namespace jian {

DistAnal::DistAnal(int cutoff, double interval, string reference_state) {
	num = NULL;
	type = NULL;
	ntLen = NULL;
	list = NULL;
	len = 0;

	for (int i = 0; i < 5; i++) {
		score[i] = 0;
	}

	this->interval = interval;
	this->cutoff = cutoff;
	this->reference_state = reference_state;
	bins = int(ceil(cutoff / interval));
	penalty = 0;

	for (int i = 0; i < 5; i++) {
		obs_parm_[i] = new int[85 * 85 * bins];
		ref_parm_[i] = new int[85 * 85 * bins];
		obs_prob_[i] = new double[85 * 85 * bins];
		ref_prob_[i] = new double[85 * 85 * bins];
		for (int j = 0; j < 85 * 85 * bins; j++) {
			obs_parm_[i][j] = 0;
			ref_parm_[i][j] = 0;
			obs_prob_[i][j] = 0;
			ref_parm_[i][j] = 0;
		}
	}
}

DistAnal::~DistAnal() {
	delete [] num;
	delete [] type;
	delete [] ntLen;
	for (int i = 0; i < len; i++) {
		delete [] list[i];
	}
	delete [] list;

	for (int i = 0; i < 5; i++) {
		delete [] obs_parm_[i];
		delete [] ref_parm_[i];
		
		delete [] obs_prob_[i];
		delete [] ref_prob_[i];
	}

	delete [] nuc_len_;
	delete [] stacking_len_;
}

void DistAnal::readRNA(const RNA &rna) {
	delete [] num;
	delete [] type;
	delete [] ntLen;
	for (int i = 0; i < len; i++) {
		delete [] list[i];
	}
	delete [] list;

	len = rna.res_nums();
	num = new int[len];
	type = new int[len];
	ntLen = new int[len];
	list = new Point *[len];

	nuc_score_ = new double[len];
	nuc_len_ = new int[len];
	for (int i = 0; i < len; i++) {
		nuc_score_[i] = 0;
		nuc_len_[i] = 0;
	}
	stacking_score_ = new double[len - 1];
	stacking_len_ = new int[len - 1];
	for (int i = 0; i < len - 1; i++) {
		stacking_score_[i] = 0;
		stacking_len_[i] = 0;
	}
	pairing_score_ = new Matr_(len, len, 0);
	pairing_len_ = new Matr_(len, len, 0);

	Point *p = new Point;
	for (int i = 0, m = 0, n = 0; i < (int) rna.chains.size(); i++) {
		for (int j = 0; j < (int) rna.chains[i].residues.size(); j++, n++) {
			m++;
			ntLen[n] = (int) rna.chains[i].residues[j].atoms.size();
			string name = rna.chains[i].residues[j].name;
			if (name == "A") {
				type[n] = 1;
			} else if (name == "U") {
				type[n] = 2;
			} else if (name == "G") {
				type[n] = 3;
			} else if (name == "C") {
				type[n] = 4;
			}
			list[n] = new Point[ntLen[n]];
			for (int k = 0; k < (int) rna.chains[i].residues[j].atoms.size(); k++) {
				list[n][k].x = rna.chains[i].residues[j].atoms[k].x;
				list[n][k].y = rna.chains[i].residues[j].atoms[k].y;
				list[n][k].z = rna.chains[i].residues[j].atoms[k].z;
				if (rna.chains[i].residues[j].atoms[0].name != "P") {
					list[n][k].type = 3;
				} else {
					list[n][k].type = 0;
				}
				if (type[n] == 1) {
					list[n][k].type += k;
				} else if (type[n] == 2) {
					list[n][k].type += 22 + k;
				} else if (type[n] == 3) {
					list[n][k].type += 42 + k;
				} else {
					list[n][k].type += 65 + k;
				}
				if (rna.chains[i].residues[j].atoms[k].name == "O5*" && n != 0) {
					double dx = list[n][k].x - p->x;
					double dy = list[n][k].y - p->y;
					double dz = list[n][k].z - p->z;
					double dist = sqrt(dx * dx + dy * dy + dz * dz);
					if (dist > 4) {
						m += 10000;
					}
				} else if (rna.chains[i].residues[j].atoms[k].name == "O3*") {
					p->x = list[n][k].x;
					p->y = list[n][k].y;
					p->z = list[n][k].z;
				}
			}
			num[n] = m;
		}
	}
}

void DistAnal::train() {
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			for (int k = 0; k < ntLen[i]; k++) {
				for (int l = 0; l < ntLen[j]; l++) {
					int type1 = list[i][k].type;
					int type2 = list[j][l].type;
					double temp = list[i][k].dist(list[j][l]);
					if (temp >= cutoff) continue;
					if (num[j] - num[i] == 0) {
						obs_parm_[0][(type1 * 85 + type2) * bins + int(temp / interval)]++;
					} else if (num[j] - num[i] == 1) {
						obs_parm_[1][(type1 * 85 + type2) * bins + int(temp / interval)]++;
					} else if (num[j] - num[i] == 2) {
						obs_parm_[2][(type1 * 85 + type2) * bins + int(temp / interval)]++;
					} else if (num[j] - num[i] == 3) {
						obs_parm_[3][(type1 * 85 + type2) * bins + int(temp / interval)]++;
					} else {
						obs_parm_[4][(type1 * 85 + type2) * bins + int(temp / interval)]++;
					}
					/*
					if (num[j] - num[i] == 1 && (((type1 >= 0 && type1 <= 11) || 
							(type1 >= 22 && type1 <= 33) || (type1 >= 42 && type1 <= 53) || 
							(type1 >= 62 && type1 <= 73)) && ((type2 >= 0 && type2 <= 11) || 
							(type2 >= 22 && type2 <= 33) || (type2 >= 42 && type2 <= 53) || 
							(type2 >= 62 && type2 <= 73)))) 
						continue;
					double temp = list[i][k].dist(&(list[j][l]));
					if (temp >= cutoff) continue;
					obsParm[(list[i][k].type * 85 + list[j][l].type) * bins + int(temp / interval)]++;
					obsParm[(list[j][l].type * 85 + list[i][k].type) * bins + int(temp / interval)]++;
					*/
					//refParm[int(temp / interval)] += 2;
				}
			}
		}
	}
}

void DistAnal::readObsParm(string filename) {
	ifstream ifile(filename.c_str());
	int temp;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85 * bins; j++) {
			ifile >> temp;
			obs_parm_[i][j] += temp;
		}
	}
	ifile.close();
	initObsProb();
	/*
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			ifile >> temp;
			if (!ifile) {
				break;
			}
			obsParm[i * bins + j] += temp;
		}
	}
	ifile.close();
	initObsProb();
	*/
}

void DistAnal::readRefParm(string filename) {
	if (reference_state == "average") {
		int *single_row_sum = new int[85 * 85];
		int *single_column_sum = new int[bins];
		for (int i = 0; i < 5; i++) {
			int total = 0;
			for (int j = 0; j < 85 * 85; j++) {
				single_row_sum[j] = 0;
				for (int k = 0; k < bins; k++) {
					if (j == 0) {
						single_column_sum[k] = 0;
					}
					single_row_sum[j] += obs_parm_[i][j * bins + k];
					single_column_sum[k] += obs_parm_[i][j * bins + k];
					total += obs_parm_[i][j * bins + k];
				}
			}
			for (int j = 0; j < 85 * 85; j++) {
				for (int k = 0; k < bins; k++) {
					ref_parm_[i][j * bins + k] = int(single_row_sum[j] * double(single_column_sum[k]) / total);
				}
			}
		}
		delete [] single_row_sum;
		delete [] single_column_sum;
		/*
		for (int i = 0; i < 85 * 85; i++) {
			single_row_sum[i] = 0;
			for (int j = 0; j < bins; j++) {
				if (i == 0) {
					single_column_sum[j] = 0;
				}
				single_row_sum[i] += obsParm[i * bins + j];
				single_column_sum[j] += obsParm[i * bins + j];
				total += obsParm[i * bins + j];
			}
		}
		for (int i = 0; i < 85 * 85; i++) {
			for (int j = 0; j < bins; j++) {
				refParm[i * bins + j] = int(single_row_sum[i] * double(single_column_sum[j]) / total);
			}
		}
		delete [] single_row_sum;
		delete [] single_column_sum;
		*/
	} else {
		ifstream ifile(filename.c_str());
		int temp;
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 85 * 85 * bins; j++) {
				ifile >> temp;
				ref_parm_[i][j] += temp;
			}
		}
		ifile.close();
	}
	initRefProb();
}

double DistAnal::scoring() {
	for (int i = 0; i < 5; i++) {
		score[i] = 0;
	}
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			for (int k = 0; k < ntLen[i]; k++) {
				for (int l = 0; l < ntLen[j]; l++) {
					int type1 = list[i][k].type;
					int type2 = list[j][l].type;
					double temp = list[i][k].dist(list[j][l]);
					if (temp >= cutoff) continue;
					double a, b, temp_score;
					if (num[j] - num[i] == 0) {
						a = obs_prob_[0][(type1 * 85 + type2) * bins + int(temp / interval)];
						b = ref_prob_[0][(type1 * 85 + type2) * bins + int(temp / interval)];
						if (a == 0 || b == 0) {
							temp_score = penalty;
						} else {
							temp_score = -log(a/b);
						}
						score[0] += temp_score;
						nuc_score_[i] += temp_score;
						nuc_len_[i]++;
					} else if (num[j] - num[i] == 1) {
						a = obs_prob_[1][(type1 * 85 + type2) * bins + int(temp / interval)];
						b = ref_prob_[1][(type1 * 85 + type2) * bins + int(temp / interval)];
						if (a == 0 || b == 0) {
							temp_score = penalty;
						} else {
							temp_score = -log(a/b);
						}
						score[1] += temp_score;
						stacking_score_[i] += temp_score;
						stacking_len_[i]++;
					} else if (num[j] - num[i] == 2) {
						a = obs_prob_[2][(type1 * 85 + type2) * bins + int(temp / interval)];
						b = ref_prob_[2][(type1 * 85 + type2) * bins + int(temp / interval)];
						if (a == 0 || b == 0) {
							temp_score = penalty;
						} else {
							temp_score = -log(a/b);
						}
						score[2] += temp_score;
					} else if (num[j] - num[i] == 3) {
						a = obs_prob_[3][(type1 * 85 + type2) * bins + int(temp / interval)];
						b = ref_prob_[3][(type1 * 85 + type2) * bins + int(temp / interval)];
						if (a == 0 || b == 0) {
							temp_score = penalty;
						} else {
							temp_score = -log(a/b);
						}
						score[3] += temp_score;
					} else {
						a = obs_prob_[4][(type1 * 85 + type2) * bins + int(temp / interval)];
						b = ref_prob_[4][(type1 * 85 + type2) * bins + int(temp / interval)];
						if (a == 0 || b == 0) {
							temp_score = penalty;
						} else {
							temp_score = -log(a/b);
						}
						score[4] += temp_score;
						pairing_score_->data[i][j] += temp_score;
						pairing_len_->data[i][j]++;
					}
				}
			}
		}
	}
//	score = score / (len * (len - 1));
	for (int i = 0; i < len; i++) {
		if (nuc_len_[i] == 0) {
			nuc_score_[i] = 0;
		} else {
			nuc_score_[i] /= nuc_len_[i];
		}
	}
	for (int i = 0; i < len - 1; i++) {
		if (stacking_len_[i] == 0) {
			stacking_score_[i] = 0;
		} else {
			stacking_score_[i] /= stacking_len_[i];
		}
	}
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (pairing_len_->data[i][j] == 0) {
				pairing_score_->data[i][j] = 0;
			} else {
				pairing_score_->data[i][j] /= pairing_len_->data[i][j];
			}
		}
	}
	return score[0] + score[1] + score[2] + score[3] + score[4];
}

void DistAnal::initObsProb() {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85; j++) {
			double line_sum = 0;
			for (int k = 0; k < bins; k++) {
				line_sum += obs_parm_[i][j * bins + k];
			}
			for (int k = 0; k < bins; k++) {
				if (line_sum == 0) {
					obs_prob_[i][j * bins + k] = 0;
				} else {
					obs_prob_[i][j * bins + k] = double(obs_parm_[i][j * bins + k]) / line_sum;
				}
			}
		}
	}
	/*
	for (int i = 0; i < 85 * 85; i++) {
		double line_sum = 0;
		for (int j = 0; j < bins; j++) {
			line_sum += obsParm[i * bins + j];
		}
		for (int j = 0; j < bins; j++) {
			obsProb[i * bins + j] = double(obsParm[i * bins + j]) / line_sum;
		}
	}
	*/
}

void DistAnal::initRefProb() {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85; j++) {
			double line_sum = 0;
			for (int k = 0; k < bins; k++) {
				line_sum += ref_parm_[i][j * bins + k];
			}
			for (int k = 0; k < bins; k++) {
				if (ref_parm_[i][j * bins + k] == 0) {
					ref_prob_[i][j * bins + k] = 0;
				} else {
					ref_prob_[i][j * bins + k] = double(ref_parm_[i][j * bins + k]) / line_sum;
				}
			}
		}
	}
	/*
	for (int i = 0; i < 85 * 85; i++) {
		double line_sum = 0;
		for (int j = 0; j < bins; j++) {
			line_sum += refParm[i * bins + j];
		}
		for (int j = 0; j < bins; j++) {
			refProb[i * bins + j] = double(refParm[i * bins + j]) / line_sum;
		}
	}
	*/
}

void DistAnal::printObsParm() {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85; j++) {
			for (int k = 0; k < bins; k++) {
				cout << obs_parm_[i][j * bins + k] << ' ';
			}
			cout << endl;
		}
	}
	/*
	int sum[bins];
	for (int i = 0; i < bins; i++) {
		sum[i] = 0;
	}
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			sum[j] += obsParm[i * bins + j];
			cout << obsParm[i * bins + j] << ' ';
		}
		cout << endl;
	}
	for (int i = 0; i < bins; i++) {
		cout << sum[i] << ' ';
	}
	cout << endl;
	*/
}

void DistAnal::printRefParm() {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85; j++) {
			for (int k = 0; k < bins; k++) {
				cout << ref_parm_[i][j * bins + k] << ' ';
			}
			cout << endl;
		}
	}
	/*
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			cout << refParm[i * bins + j] << ' ';
		}
		cout << endl;
	}
	*/
}

void DistAnal::printObsProb() {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85; j++) {
			for (int k = 0; k < bins; k++) {
				cout << obs_prob_[i][j * bins + k] << ' ';
			}
			cout << endl;
		}
	}
	/*
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			cout << obsProb[i * bins + j] << ' ';
		}
		cout << endl;
	}
	*/
}

void DistAnal::printRefProb() {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 85 * 85; j++) {
			for (int k = 0; k < bins; k++) {
				cout << ref_prob_[i][j * bins + k] << ' ';
			}
			cout << endl;
		}
	}
	/*
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			cout << refProb[i * bins + j] << ' ';
		}
		cout << endl;
	}
	*/
}

double *DistAnal::getScore() {
	return score;
}

} /// namespace jian


