#ifndef JIAN_NUC2D_SEQ2QUA_H
#define JIAN_NUC2D_SEQ2QUA_H

#include <util/util.h>

namespace jian {
namespace nuc2d {

namespace seq2qua {
struct MyHash {
    std::size_t operator ()(const std::vector<int> &s) const {
        long hash = 5381;  
        for(int i = 0; i < s.size(); i++)  {  
            hash = ((hash << 5) + hash) + s[i];  
        }  
        return hash;
    }
};
}

class Seq2Qua {
public:
    typedef std::tuple<int, int, int, int> Pair;
    typedef std::vector<Pair> Pairs;
    typedef std::vector<std::pair<Pairs, double>> PairLists;

    std::unordered_map<std::vector<int>, double, seq2qua::MyHash> _scores;

    std::string _seq;
    std::map<char, int> _convert{{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
    int _min_hairpin_size = 4;
    double _cutoff = 1;

//    void operator ()(std::string seq);
//    double get(const std::vector<int> &sequence);
//    double score(int, int, int, int);
//    PairLists backtrack(const std::vector<int> &sequence, double cutoff);

    void operator ()(std::string seq) {
        _seq = seq;
        int len = seq.size();
        std::vector<int> vec(len);
        std::iota(vec.begin(), vec.end(), 0);
        get(vec);

        auto pair_lists = backtrack(vec, _cutoff);
        for (auto &&pair_list: pair_lists) {
            for (auto &&pair: pair_list.first) {
                std::cout << std::get<0>(pair) << '-' << std::get<1>(pair) << '-' << std::get<2>(pair) << '-' << std::get<3>(pair) << ' ';
            }
            std::cout << pair_list.second << std::endl;
        }
    }

    double get(const std::vector<int> &seq) {
        if (seq.size() < 4) {
            return 0;
        }

        if (_scores.count(seq)) {
            return _scores[seq];
        }

        std::vector<int> temp_vec;
        std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_vec));
        double min_score = get(temp_vec);

        for (int i = 0; i < seq.size() - 3; i++) {
            for (int j = i + 1; j < seq.size() - 2; j++) {
                for (int k = j + 1; k < seq.size() - 1; k++) {
                    if (seq[j] - seq[i] < _min_hairpin_size + 1 || seq[k] - seq[j] < _min_hairpin_size + 1 || seq.back() - seq[k] < _min_hairpin_size + 1) {
                        continue;
                    }
                    std::vector<int> temp_seq1, temp_seq2;
                    std::copy(seq.begin(), std::next(seq.begin(), i), std::back_inserter(temp_seq1));
                    std::copy(std::next(seq.begin(), j + 1), std::next(seq.begin(), k), std::back_inserter(temp_seq1));
                    std::copy(std::next(seq.begin(), i + 1), std::next(seq.begin(), j), std::back_inserter(temp_seq2));
                    std::copy(std::next(seq.begin(), k + 1), std::prev(seq.end(), 1), std::back_inserter(temp_seq2));
                    double temp = get(temp_seq1) + get(temp_seq2) + score(seq[i], seq[j], seq[k], seq.back());
                    if (temp < min_score) {
                        min_score = temp;
                    }
                }
            }
        }

        _scores[seq] = min_score;
        return min_score;
    }

    double score(int a, int b, int c, int d) {
        int m = _convert[_seq[a]] * _convert[_seq[b]] * _convert[_seq[c]] * _convert[_seq[d]];
        if (m == 16) {
            return -1;
        } else {
            return 0;
        }
    }

    PairLists backtrack(const std::vector<int> &seq, double cutoff) {
        PairLists vec;

        if (seq.size() < 4) {
            return vec;
        }

        double temp_score = get(seq);

        std::vector<int> temp_seq;
        std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_seq));
        double s = get(temp_seq);
        if (s - temp_score <= cutoff) {
            auto temp_vec = backtrack(temp_seq, cutoff - s + temp_score);
            std::copy(temp_vec.begin(), temp_vec.end(), std::back_inserter(vec));
        }

        for (int i = 0; i < seq.size() - 3; i++) {
            for (int j = i + 1; j < seq.size() - 2; j++) {
                for (int k = j + 1; k < seq.size() - 1; k++) {
                    if (seq[j] - seq[i] < _min_hairpin_size + 1 || seq[k] - seq[j] < _min_hairpin_size + 1 || seq.back() - seq[k] < _min_hairpin_size + 1) {
                        continue;
                    }
                    std::vector<int> temp_seq1, temp_seq2;
                    std::copy(seq.begin(), std::next(seq.begin(), i), std::back_inserter(temp_seq1));
                    std::copy(std::next(seq.begin(), j + 1), std::next(seq.begin(), k), std::back_inserter(temp_seq1));
                    std::copy(std::next(seq.begin(), i + 1), std::next(seq.begin(), j), std::back_inserter(temp_seq2));
                    std::copy(std::next(seq.begin(), k + 1), std::prev(seq.end(), 1), std::back_inserter(temp_seq2));
                    auto en = score(seq[i], seq[j], seq[k], seq.back());
                    double temp = get(temp_seq1) + get(temp_seq2) + en;
                    if (temp == temp_score && en == -1) {
                        double new_cutoff = cutoff - temp + temp_score;
                        auto pair = std::make_tuple(seq[i], seq[j], seq[k], seq.back());
                        Pairs pair_list;
                        pair_list.push_back(pair);
                        auto pair_lists1 = backtrack(temp_seq1, new_cutoff);
                        auto pair_lists2 = backtrack(temp_seq2, new_cutoff);
                        if (pair_lists1.size() != 0 && pair_lists2.size() != 0) {
                            for (int ii = 0; ii < pair_lists1.size(); ii++) {
                                for (int jj = 0; jj < pair_lists2.size(); jj++) {
                                    double s = pair_lists1[ii].second + pair_lists2[jj].second + en;
                                    if (s - temp_score <= new_cutoff) {
                                        Pairs temp_pairs = pair_list;
                                        std::copy(pair_lists1[ii].first.begin(), pair_lists1[ii].first.end(), std::back_inserter(temp_pairs));
                                        std::copy(pair_lists2[jj].first.begin(), pair_lists2[jj].first.end(), std::back_inserter(temp_pairs));
                                        vec.push_back(std::make_pair(temp_pairs, s));
                                    }
                                }
                            }
                        } else if (pair_lists1.size() != 0) {
                            for (int ii = 0; ii < pair_lists1.size(); ii++) {
                                double s = pair_lists1[ii].second + en;
                                if (s - temp_score <= new_cutoff) {
                                    Pairs temp_pairs = pair_list;
                                    std::copy(pair_lists1[ii].first.begin(), pair_lists1[ii].first.end(), std::back_inserter(temp_pairs));
                                    vec.push_back(std::make_pair(temp_pairs, s));
                                }
                            }
                        } else if (pair_lists2.size() != 0) {
                            for (int jj = 0; jj < pair_lists2.size(); jj++) {
                                double s = pair_lists2[jj].second + en;
                                if (s - temp_score <= new_cutoff) {
                                    Pairs temp_pairs = pair_list;
                                    std::copy(pair_lists2[jj].first.begin(), pair_lists2[jj].first.end(), std::back_inserter(temp_pairs));
                                    vec.push_back(std::make_pair(temp_pairs, s));
                                }
                            }
                        } else {
                            if (en - temp_score <= new_cutoff) {
                                vec.push_back(std::make_pair(pair_list, en));
                            }
                        }
                    }
                }
            }
        }
        return vec;
    }

};

} /// namespace nuc2d
} /// namespace jian

#endif
