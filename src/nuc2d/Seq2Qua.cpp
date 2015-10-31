#include "Seq2Qua.h"

namespace jian {
namespace nuc2d {

std::vector<std::vector<std::tuple<int, int, int, int>>> Seq2Qua::operator ()(std::string seq) {
    _seq = seq;
    int len = seq.size();
    std::vector<int> vec(len);
    std::iota(vec.begin(), vec.end(), 0);
    get(vec);

    auto pair_lists = backtrack(vec);
    for (auto &&pair_list: pair_lists) {
        for (auto &&pair: pair_list) {
            std::cout << std::get<0>(pair) << '-' << std::get<1>(pair) << '-' << std::get<2>(pair) << '-' << std::get<3>(pair) << ' ';
        }
        std::cout << std::endl;
    }
    return pair_lists;
}

double Seq2Qua::get(const std::vector<int> &seq) {
    if (seq.size() < 4) {
        return 0;
    }

    if (_scores.count(seq)) {
        return _scores[seq];
    }

    std::vector<int> temp_vec;
    std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_vec));
    double max_score = get(temp_vec);

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
                if (temp > max_score) {
                    max_score = temp;
                }
            }
        }
    }

    _scores[seq] = max_score;
    return max_score;
}

double Seq2Qua::score(int a, int b, int c, int d) {
    int m = _convert[_seq[a]] * _convert[_seq[b]] * _convert[_seq[c]] * _convert[_seq[d]];
    if (m == 16) {
        return 1;
    } else {
        return 0;
    }
}

std::vector<std::vector<std::tuple<int, int, int, int>>> Seq2Qua::backtrack(const std::vector<int> &seq) {
    std::vector<std::vector<std::tuple<int, int, int, int>>> vec;

    if (seq.size() < 4) {
        return vec;
    }

    double temp_score = get(seq);

    std::vector<int> temp_seq;
    std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_seq));
    if (temp_score == get(temp_seq)) {
        vec = backtrack(temp_seq);
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
                if (temp == temp_score && en == 1) {
                    auto pair = std::make_tuple(seq[i], seq[j], seq[k], seq.back());
                    std::vector<std::tuple<int, int, int, int>> pair_list;
                    pair_list.push_back(pair);
                    std::vector<std::vector<std::tuple<int, int, int, int>>> pair_lists;
                    pair_lists.push_back(pair_list);
                    auto pair_lists1 = backtrack(temp_seq1);
                    auto pair_lists2 = backtrack(temp_seq2);
                    if (pair_lists1.size() != 0) {
                        std::vector<std::vector<std::tuple<int, int, int, int>>> temp_lists;
                        for (int i = 0; i < pair_lists.size(); i++) {
                            for (int j = 0; j < pair_lists1.size(); j++) {
                                auto temp_pair_list = pair_lists[i];
                                std::copy(pair_lists1[j].begin(), pair_lists1[j].end(), std::back_inserter(temp_pair_list));
                                temp_lists.push_back(temp_pair_list);
                            }
                        }
                        std::swap(pair_lists, temp_lists);
                    }
                    if (pair_lists2.size() != 0) {
                        std::vector<std::vector<std::tuple<int, int, int, int>>> temp_lists;
                        for (int i = 0; i < pair_lists.size(); i++) {
                            for (int j = 0; j < pair_lists2.size(); j++) {
                                auto temp_pair_list = pair_lists[i];
                                std::copy(pair_lists2[j].begin(), pair_lists2[j].end(), std::back_inserter(temp_pair_list));
                                temp_lists.push_back(temp_pair_list);
                            }
                        }
                        std::swap(pair_lists, temp_lists);
                    }
                    std::copy(pair_lists.begin(), pair_lists.end(), std::back_inserter(vec));
                }
            }
        }
    }
    return vec;
}

} /// namespace nuc2d
} /// namespace jian
