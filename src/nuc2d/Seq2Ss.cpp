#include "Seq2Ss.h"

namespace jian {
namespace nuc2d {

void Seq2Ss::operator ()(std::string seq) {
    _len = seq.size();
    _seq = seq;
    std::map<char, int> temp_map{{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
    std::transform(_seq.begin(), _seq.end(), std::back_inserter(_types), [&](char c){return temp_map[c];});
    _mms = MatrixXf::Zero(_len, _len);
    for (int i = 0; i < _len - 1 - _min_hairpin_size; i++) {
        for (int j = 0; j < _len - 1 - _min_hairpin_size - i; j++) {
            _mms(j, i + j + 1 + _min_hairpin_size) = mms(j, i + j + 1 + _min_hairpin_size);
        }
    }
    auto all_pairs = backtrack(0, _len - 1, _cutoff);
    for (auto &&pairs: all_pairs) {
        std::string s(_seq.size(), '.');
        for (auto &&pair: pairs.first) {
            s[pair.first] = '(';
            s[pair.second] = ')';
        }
        std::cout << s << ' ' << pairs.second << std::endl;
    }
}

double Seq2Ss::mms(int m, int n) {
    if (n <= m + _min_hairpin_size) {
        return 0;
    }
    double min_score = get_mms(m, n - 1);
    for (int i = m; i < n - _min_hairpin_size; i++) {
        double temp = get_mms(m, i - 1) + get_mms(i + 1, n - 1) + score(i, n);
        if (temp < min_score) {
            min_score = temp;
        }
    }
    return min_score;
}

double Seq2Ss::score(int m, int n) {
    int a = _types[m] + _types[n];
    int b = _types[m] * _types[n];
    if (a == 1 || a == 5 || b == 2) {
        return -1;
    } else {
        return 0;
    }
}

Seq2Ss::PairLists Seq2Ss::backtrack(int m, int n, double cutoff) {
    PairLists temp_vec;
    if (n <= m + _min_hairpin_size) {
        return temp_vec;;
    }
//    double min_score = get_mms(m, n);
//    if (-min_score <= cutoff) {
//        temp_vec.push_back(std::make_pair(Pairs(), 0));
//    }
    double temp_score = get_mms(m, n - 1) - min_score;
    if (temp_score <= cutoff) {
        auto vec = backtrack(m, n - 1, cutoff - temp_score);
        std::copy(vec.begin(), vec.end(), std::back_inserter(temp_vec));
    }
    for (int i = m; i < n - _min_hairpin_size; i++) {
        temp_score = get_mms(m, i - 1) + get_mms(i + 1, n - 1) + score(i, n) - min_score;
        if (temp_score <= cutoff) {
            auto pairs1 = backtrack(m, i - 1, cutoff - temp_score);
            auto pairs2 = backtrack(i + 1, n - 1, cutoff - temp_score);
            Pairs pairs;
            if (score(i, n) == -1) {
                pairs.push_back(make_pair(i, n));
            }
            if (pairs1.size() != 0 && pairs2.size() != 0) {
                for (int j = 0; j < pairs1.size(); j++) {
                    for (int k = 0; k < pairs2.size(); k++) {
                        double temp = pairs1[j].second + pairs2[k].second + score(i, n);
                        if (temp - min_score <= cutoff) {
                            auto temp_pairs = pairs;
                            std::copy(pairs1[j].first.begin(), pairs1[j].first.end(), std::back_inserter(temp_pairs));
                            std::copy(pairs2[k].first.begin(), pairs2[k].first.end(), std::back_inserter(temp_pairs));
                            temp_vec.push_back(std::make_pair(temp_pairs, temp));
                        }
                    }
                }
            } else if (pairs1.size() != 0) {
                for (int j = 0; j < pairs1.size(); j++) {
                    double temp = pairs1[j].second + score(i, n);
                    if (temp - min_score <= cutoff) {
                        auto temp_pairs = pairs;
                        std::copy(pairs1[j].first.begin(), pairs1[j].first.end(), std::back_inserter(temp_pairs));
                        temp_vec.push_back(std::make_pair(temp_pairs, temp));
                    }
                }
            } else if (pairs2.size() != 0) {
                for (int j = 0; j < pairs2.size(); j++) {
                    double temp = pairs2[j].second + score(i, n);
                    if (temp - min_score <= cutoff) {
                        auto temp_pairs = pairs;
                        std::copy(pairs2[j].first.begin(), pairs2[j].first.end(), std::back_inserter(temp_pairs));
                        temp_vec.push_back(std::make_pair(temp_pairs, temp));
                    }
                }
            } else {
                temp_vec.push_back(std::make_pair(pairs, score(i, n)));
            }
        }
    }
    return temp_vec;
}

double Seq2Ss::get_mms(int m, int n) {
    if (n <= m + _min_hairpin_size) {
        return 0;
    } else {
        return _mms(m, n);
    }
    
}
        
} /// namespace nuc2d
} /// namespace jian
