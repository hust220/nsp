#include "Seq2Tri.h"

namespace jian {
namespace nuc2d {

void Seq2Tri::operator ()(std::string seq) {
    _seq = seq;
    int len = seq.size();
    std::vector<int> vec(len);
    std::iota(vec.begin(), vec.end(), 0);
    get(vec);

    auto pair_lists = backtrack(vec, _cutoff);
    for (auto &&pair_list: pair_lists) {
        for (auto &&pair: pair_list.first) {
            std::cout << std::get<0>(pair) << '-' << std::get<1>(pair) << '-' << std::get<2>(pair) << ' ';
        }
        std::cout << pair_list.second << std::endl;
    }
}

double Seq2Tri::get(const std::vector<int> &seq) {
    if (seq.size() < 3) {
        return 0;
    }

    if (_scores.count(seq)) {
        return _scores[seq];
    }

    std::vector<int> temp_vec;
    std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_vec));
    double min_score = get(temp_vec);

    for (int i = 0; i < seq.size() - 2; i++) {
        for (int j = i + 1; j < seq.size() - 1; j++) {
            if (seq[j] - seq[i] < _min_hairpin_size + 1 || seq.back() - seq[j] < _min_hairpin_size + 1) {
                continue;
            }
            std::vector<int> temp_seq1, temp_seq2;
            std::copy(seq.begin(), std::next(seq.begin(), i), std::back_inserter(temp_seq1));
            std::copy(std::next(seq.begin(), j + 1), std::prev(seq.end(), 1), std::back_inserter(temp_seq1));
            std::copy(std::next(seq.begin(), i + 1), std::next(seq.begin(), j), std::back_inserter(temp_seq2));
            double temp = get(temp_seq1) + get(temp_seq2) + score(seq[i], seq[j], seq.back());
            if (temp < min_score) {
                min_score = temp;
            }
        }
    }

    _scores[seq] = min_score;
    return min_score;
}

double Seq2Tri::score(int a, int b, int c) {
    int m = _convert[_seq[a]] * _convert[_seq[b]] * _convert[_seq[c]];
    if (m == 12 || m == 18) {
        return -1;
    } else {
        return 0;
    }
}

Seq2Tri::PairLists Seq2Tri::backtrack(const std::vector<int> &seq, double cutoff) {
    PairLists vec;

    if (seq.size() < 3) {
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

    for (int i = 0; i < seq.size() - 2; i++) {
        for (int j = i + 1; j < seq.size() - 1; j++) {
            if (seq[j] - seq[i] < _min_hairpin_size + 1 || seq.back() - seq[j] < _min_hairpin_size + 1) {
                continue;
            }
            std::vector<int> temp_seq1, temp_seq2;
            std::copy(seq.begin(), std::next(seq.begin(), i), std::back_inserter(temp_seq1));
            std::copy(std::next(seq.begin(), j + 1), std::prev(seq.end(), 1), std::back_inserter(temp_seq1));
            std::copy(std::next(seq.begin(), i + 1), std::next(seq.begin(), j), std::back_inserter(temp_seq2));
            auto en = score(seq[i], seq[j], seq.back());
            double temp = get(temp_seq1) + get(temp_seq2) + en;
            if (temp - temp_score <= cutoff && en == -1) {
                auto pair = std::make_tuple(seq[i], seq[j], seq.back());
                Pairs pair_list;
                pair_list.push_back(pair);
                double new_cutoff = cutoff - (temp - temp_score);
                auto pair_lists1 = backtrack(temp_seq1, new_cutoff);
                auto pair_lists2 = backtrack(temp_seq2, new_cutoff);
                if (pair_lists1.size() != 0 && pair_lists2.size() != 0) {
                    for (int k = 0; k < pair_lists1.size(); k++) {
                        for (int l = 0; l < pair_lists2.size(); l++) {
                            double s = pair_lists1[k].second + pair_lists2[l].second + en;
                            if (s - temp_score <= new_cutoff) {
                                Pairs temp_pairs = pair_list;
                                std::copy(pair_lists1[k].first.begin(), pair_lists1[k].first.end(), std::back_inserter(temp_pairs));
                                std::copy(pair_lists2[l].first.begin(), pair_lists2[l].first.end(), std::back_inserter(temp_pairs));
                                vec.push_back(std::make_pair(temp_pairs, s));
                            }
                        }
                    }
                } else if (pair_lists1.size() != 0) {
                    for (int k = 0; k < pair_lists1.size(); k++) {
                        double s = pair_lists1[k].second + en;
                        if (s - temp_score <= new_cutoff) {
                            Pairs temp_pairs = pair_list;
                            std::copy(pair_lists1[k].first.begin(), pair_lists1[k].first.end(), std::back_inserter(temp_pairs));
                            vec.push_back(std::make_pair(temp_pairs, s));
                        }
                    }
                } else if (pair_lists2.size() != 0) {
                    for (int l = 0; l < pair_lists2.size(); l++) {
                        double s = pair_lists2[l].second + en;
                        if (s - temp_score <= new_cutoff) {
                            Pairs temp_pairs = pair_list;
                            std::copy(pair_lists2[l].first.begin(), pair_lists2[l].first.end(), std::back_inserter(temp_pairs));
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
    return vec;
}

} /// namespace nuc2d
} /// namespace jian
