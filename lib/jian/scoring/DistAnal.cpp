#include "DistAnal.h"

namespace jian {
namespace scoring {

DistAnal::DistAnal(int cutoff, double interval, std::string ref_state) {
    _interval = interval;
    _cutoff = cutoff;
    _ref_state = ref_state;
    _bins = int(ceil(_cutoff / _interval));

    _obs_parm = MatrixXd::Zero(5, 85 * 85 * _bins);
    _ref_parm = MatrixXd::Zero(5, 85 * 85 * _bins);
    _obs_prob = MatrixXf::Zero(5, 85 * 85 * _bins);
    _ref_prob = MatrixXf::Zero(5, 85 * 85 * _bins);
}

void DistAnal::read_mol(const Model &rna) {
    _len = rna.res_nums();
    _num.resize(_len);
    _type.resize(_len);
    _ntLen.resize(_len);
    _list.resize(_len);
    _nuc_score = VectorXf::Zero(_len);
    _nuc_len = VectorXd::Zero(_len);
    _stacking_score = VectorXf::Zero(_len - 1);
    _stacking_len = VectorXd::Zero(_len - 1);
    _pairing_score = MatrixXf::Zero(_len, _len);
    _pairing_len = MatrixXd::Zero(_len, _len);

    Point p;
    std::map<std::string, int> name_type_map{{"A", 1}, {"U", 2}, {"G", 3}, {"C", 4}};
    for (int i = 0, m = 0, n = 0; i < (int) rna.chains.size(); i++) {
        for (int j = 0; j < (int) rna.chains[i].residues.size(); j++, n++) {
            m++;
            _ntLen[n] = (int) rna.chains[i].residues[j].atoms.size();
            string name = rna.chains[i].residues[j].name;
            _type[n] = name_type_map[name];
            _list[n].resize(_ntLen[n]);
            for (int k = 0; k < (int) rna.chains[i].residues[j].atoms.size(); k++) {
                _list[n][k].x = rna.chains[i].residues[j].atoms[k].x;
                _list[n][k].y = rna.chains[i].residues[j].atoms[k].y;
                _list[n][k].z = rna.chains[i].residues[j].atoms[k].z;
                if (rna.chains[i].residues[j].atoms[0].name != "P") {
                    _list[n][k].type = 3;
                } else {
                    _list[n][k].type = 0;
                }
                if (_type[n] == 1) {
                    _list[n][k].type += k;
                } else if (_type[n] == 2) {
                    _list[n][k].type += 22 + k;
                } else if (_type[n] == 3) {
                    _list[n][k].type += 42 + k;
                } else {
                    _list[n][k].type += 65 + k;
                }
                if (rna.chains[i].residues[j].atoms[k].name == "O5*" && n != 0) {
                    if (_list[n][k].dist(p) > 4) {
                        m += 10000;
                    }
                } else if (rna.chains[i].residues[j].atoms[k].name == "O3*") {
                    p = _list[n][k];
                }
            }
            _num[n] = m;
        }
    }
}

void DistAnal::train() {
    for (int i = 0; i < _len; i++) {
        for (int j = i; j < _len; j++) {
            for (int k = 0; k < _ntLen[i]; k++) {
                for (int l = 0; l < _ntLen[j]; l++) {
                    int type1 = _list[i][k].type;
                    int type2 = _list[j][l].type;
                    double temp = _list[i][k].dist(_list[j][l]);
                    if (temp >= _cutoff) continue;
                    if (_num[j] - _num[i] == 0) {
                        _obs_parm(0, (type1 * 85 + type2) * _bins + int(temp / _interval))++;
                    } else if (_num[j] - _num[i] == 1) {
                        _obs_parm(1, (type1 * 85 + type2) * _bins + int(temp / _interval))++;
                    } else if (_num[j] - _num[i] == 2) {
                        _obs_parm(2, (type1 * 85 + type2) * _bins + int(temp / _interval))++;
                    } else if (_num[j] - _num[i] == 3) {
                        _obs_parm(3, (type1 * 85 + type2) * _bins + int(temp / _interval))++;
                    } else {
                        _obs_parm(4, (type1 * 85 + type2) * _bins + int(temp / _interval))++;
                    }
                }
            }
        }
    }
}

void DistAnal::read_obs_parm(string filename) {
    ifstream ifile(filename.c_str());
    int temp;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 85 * 85 * _bins; j++) {
            ifile >> temp;
            _obs_parm(i, j) += temp;
        }
    }
    ifile.close();
    init_obs_prob();
}

void DistAnal::read_ref_parm(string filename) {
    if (_ref_state == "average") {
        int *single_row_sum = new int[85 * 85];
        int *single_column_sum = new int[_bins];
        for (int i = 0; i < 5; i++) {
            int total = 0;
            for (int j = 0; j < 85 * 85; j++) {
                single_row_sum[j] = 0;
                for (int k = 0; k < _bins; k++) {
                    if (j == 0) {
                        single_column_sum[k] = 0;
                    }
                    single_row_sum[j] += _obs_parm(i, j * _bins + k);
                    single_column_sum[k] += _obs_parm(i, j * _bins + k);
                    total += _obs_parm(i, j * _bins + k);
                }
            }
            for (int j = 0; j < 85 * 85; j++) {
                for (int k = 0; k < _bins; k++) {
                    _ref_parm(i, j * _bins + k) = int(single_row_sum[j] * double(single_column_sum[k]) / total);
                }
            }
        }
        delete [] single_row_sum;
        delete [] single_column_sum;
    } else {
        ifstream ifile(filename.c_str());
        int temp;
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 85 * 85 * _bins; j++) {
                ifile >> temp;
                _ref_parm(i, j) += temp;
            }
        }
        ifile.close();
    }
    init_ref_prob();
}

double DistAnal::operator ()(const Model &model) {
    read_mol(model);
    _scores = VectorXf::Zero(5);
    for (int i = 0; i < _len; i++) {
        for (int j = i; j < _len; j++) {
            for (int k = 0; k < _ntLen[i]; k++) {
                for (int l = 0; l < _ntLen[j]; l++) {
                    int type1 = _list[i][k].type;
                    int type2 = _list[j][l].type;
                    double temp = _list[i][k].dist(_list[j][l]);
                    if (temp >= _cutoff) continue;
                    double a, b, temp_score;
                    if (_num[j] - _num[i] == 0) {
                        a = _obs_prob(0, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        b = _ref_prob(0, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        if (a == 0 || b == 0) {
                            temp_score = _penalty;
                        } else {
                            temp_score = -log(a/b);
                        }
                        _scores[0] += temp_score;
                        _nuc_score[i] += temp_score;
                        _nuc_len[i]++;
                    } else if (_num[j] - _num[i] == 1) {
                        a = _obs_prob(1, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        b = _ref_prob(1, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        if (a == 0 || b == 0) {
                            temp_score = _penalty;
                        } else {
                            temp_score = -log(a/b);
                        }
                        _scores[1] += temp_score;
                        _stacking_score[i] += temp_score;
                        _stacking_len[i]++;
                    } else if (_num[j] - _num[i] == 2) {
                        a = _obs_prob(2, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        b = _ref_prob(2, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        if (a == 0 || b == 0) {
                            temp_score = _penalty;
                        } else {
                            temp_score = -log(a/b);
                        }
                        _scores[2] += temp_score;
                    } else if (_num[j] - _num[i] == 3) {
                        a = _obs_prob(3, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        b = _ref_prob(3, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        if (a == 0 || b == 0) {
                            temp_score = _penalty;
                        } else {
                            temp_score = -log(a/b);
                        }
                        _scores[3] += temp_score;
                    } else {
                        a = _obs_prob(4, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        b = _ref_prob(4, (type1 * 85 + type2) * _bins + int(temp / _interval));
                        if (a == 0 || b == 0) {
                            temp_score = _penalty;
                        } else {
                            temp_score = -log(a/b);
                        }
                        _scores[4] += temp_score;
                        _pairing_score(i, j) += temp_score;
                        _pairing_len(i, j)++;
                    }
                }
            }
        }
    }
    for (int i = 0; i < _len; i++) {
        if (_nuc_len[i] == 0) {
            _nuc_score[i] = 0;
        } else {
            _nuc_score[i] /= _nuc_len[i];
        }
    }
    for (int i = 0; i < _len - 1; i++) {
        if (_stacking_len[i] == 0) {
            _stacking_score[i] = 0;
        } else {
            _stacking_score[i] /= _stacking_len[i];
        }
    }
    for (int i = 0; i < _len; i++) {
        for (int j = i; j < _len; j++) {
            if (_pairing_len(i, j) == 0) {
                _pairing_score(i, j) = 0;
            } else {
                _pairing_score(i, j) /= _pairing_len(i, j);
            }
        }
    }
    return _scores.sum();
}

void DistAnal::init_obs_prob() {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 85 * 85; j++) {
            double line_sum = 0;
            for (int k = 0; k < _bins; k++) line_sum += _obs_parm(i, j * _bins + k);
            for (int k = 0; k < _bins; k++) {
                if (line_sum == 0) {
                    _obs_prob(i, j * _bins + k) = 0;
                } else {
                    _obs_prob(i, j * _bins + k) = double(_obs_parm(i, j * _bins + k)) / line_sum;
                }
            }
        }
    }
}

void DistAnal::init_ref_prob() {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 85 * 85; j++) {
            double line_sum = 0;
            for (int k = 0; k < _bins; k++) line_sum += _ref_parm(i, j * _bins + k);
            for (int k = 0; k < _bins; k++) {
                if (_ref_parm(i, j * _bins + k) == 0) _ref_prob(i, j * _bins + k) = 0;
                else _ref_prob(i, j * _bins + k) = double(_ref_parm(i, j * _bins + k)) / line_sum;
            }
        }
    }
}

double DistAnal::score() {
    return _scores.sum();
}

} /// namespace scoring
} /// namespace jian


