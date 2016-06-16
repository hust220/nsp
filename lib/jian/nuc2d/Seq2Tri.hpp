#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>
#include <string>
#include <unordered_map>
#include "../utils/ls.hpp"

namespace jian {
namespace nuc2d {

namespace seq2tri {
struct MyHash {
    std::size_t operator ()(const std::vector<int> &s) const {
        long hash = 5381;  
        for(int i = 0; i < s.size(); i++) hash = ((hash << 5) + hash) + s[i];
        return hash;
    }
};
}

class Seq2Tri {
public:
    typedef std::tuple<int, int, int> Tuple;
    typedef std::vector<Tuple> Tuples;
    typedef std::vector<std::vector<int>> TupleList;
    typedef std::pair<TupleList, double> TupleInfo;
    typedef std::vector<TupleInfo> InfoList;

    std::unordered_map<std::vector<int>, InfoList, seq2tri::MyHash> _info;

    std::string _seq;
    std::vector<int> _types;
    std::vector<std::tuple<int, int, int>> _pairs;
    std::map<char, int> _convert{{'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3}};
    std::hash<std::string> hash_fn;
    int _min_hairpin_size;
    double _cutoff;
    Eigen::Matrix4f _pair_energy;
    Eigen::Matrix4f _stack_energy;

//    Seq2Tri();
//    InfoList operator ()(std::string seq);
//    std::string dbn(const TupleInfo &info, int len);
//    InfoList sub_info(const std::vector<int> &sequence, double cutoff);
//    InfoList best_info(const std::vector<int> &sequence);
//    TupleInfo strand_info(const std::vector<int>  &sequence);
//    TupleInfo score(const TupleInfo &info1, const TupleInfo &info2, int a, int b, int c);
//    TupleInfo score(const TupleInfo &info, int m);
//    std::pair<double, bool> tuple_energy(char i, char j, char k);
//    std::pair<double, bool> pair_energy(char m, char n);
//    double stack_energy(char i, char j);
//
    Seq2Tri() {
        _cutoff = 0.3;
        _min_hairpin_size = 4;
        _pair_energy <<  -45.94, -484.06,  -45.94,   -45.94,
                        -491.53,  -45.94,    -200,   -45.94,
                         -45.94,    -200,  -45.94,  -679.05,
                         -45.94,  -45.94, -663.93,   -45.94;
        _stack_energy << -625.57, -578.22, -772.93, -650.16,
                         -700.71, -580.50, -753.00, -525.07,
                         -742.75, -700.00, -837.22, -810.87,
                         -642.19, -735.33, -851.11, -509.75;
    }

    InfoList operator ()(std::string seq) {
        _seq = seq;
        int len = seq.size();
    //    std::vector<int> vec(len); std::iota(vec.begin(), vec.end(), 0);
        auto vec = range<std::vector<int>>(0, len);
        auto best_info_list = best_info(vec);

        auto sub_info_list = sub_info(vec, -_cutoff * best_info_list[0].second);
        std::sort(sub_info_list.begin(), sub_info_list.end(), [](const TupleInfo &info1, const TupleInfo &info2){
            return info1.second < info2.second;});
        return sub_info_list;
    }

    std::string dbn(const TupleInfo &info, int len) {
        std::string ss(len, '.');
        for (int i = 0; i < info.first.size(); i++) {
            ss[info.first[i][0]] = '<'; ss[info.first[i][1]] = 'x'; ss[info.first[i][2]] = '>';
        }
        return ss;
    }

    InfoList best_info(const std::vector<int> &seq) {
        InfoList info_list;

        if (seq.size() < 3) {
            info_list.push_back(strand_info(seq));
            return info_list;
        }

        if (_info.count(seq)) {
            return _info[seq];
        }

        /// add all the best information of first n - 1 nucleotides
        InfoList temp_info_list;
        std::vector<int> temp_vec;
        std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_vec));
        auto temp_best_info = best_info(temp_vec);
        for (auto &&info: temp_best_info) {
            temp_info_list.push_back(score(info, seq.back()));
        }

        for (int i = 0; i < seq.size() - 2; i++) {
            for (int j = i + 1; j < seq.size() - 1; j++) {
                if (seq[j] - seq[i] < _min_hairpin_size + 1 || seq.back() - seq[j] < _min_hairpin_size + 1 ||
                        tuple_energy(_seq[seq[i]], _seq[seq[j]], _seq[seq.back()]).second == false) {
                    continue;
                }
                /// construct fragment 1 and fragment 2
                std::vector<int> temp_seq1, temp_seq2;
                std::copy(seq.begin(), std::next(seq.begin(), i), std::back_inserter(temp_seq1));
                std::copy(std::next(seq.begin(), j + 1), std::prev(seq.end(), 1), std::back_inserter(temp_seq1));
                std::copy(std::next(seq.begin(), i + 1), std::next(seq.begin(), j), std::back_inserter(temp_seq2));

                /// calculate best information
                auto best_info1 = best_info(temp_seq1);
                auto best_info2 = best_info(temp_seq2);
                for (auto &&info1: best_info1) {
                    for (auto &&info2: best_info2) {
                        temp_info_list.push_back(score(info1, info2, seq[i], seq[j], seq.back()));
                    }
                }
            }
        }

        std::vector<int> score_list;
        score_list.push_back(0);
        for (int i = 1; i < temp_info_list.size(); i++) {
            if (temp_info_list[i].first.size() == 0) continue;
            double en1 = temp_info_list[i].second;
            double en2 = temp_info_list[score_list[0]].second;
            if (en1 < en2) {
                score_list.clear();
                score_list.push_back(i);
            } else if (en1 == en2) {
                score_list.push_back(i);
            }
        }
        for (auto &&i: score_list) {
            info_list.push_back(temp_info_list[i]);
        }

        _info[seq] = info_list;
        return info_list;
    }

    TupleInfo strand_info(const std::vector<int> &seq) {
        TupleInfo tuple_info;
        int len = seq.size();
        tuple_info.first.resize(len);
        tuple_info.second = 0;
        return tuple_info;
    }

    TupleInfo score(const TupleInfo &info, int m) {
        int len = info.first.size() + 1;
        TupleInfo tuple_info;
        tuple_info.first.resize(len);
        std::copy(info.first.begin(), info.first.end(), tuple_info.first.begin());
        tuple_info.second = info.second + 0;
        return tuple_info;
    }

    TupleInfo score(const TupleInfo &info1, const TupleInfo &info2, int a, int b, int c) {
        int len1 = info1.first.size() - c + b + 1;
        int len2 = info2.first.size();
        int len3 = c - b - 1;

        TupleInfo tuple_info;
        tuple_info.first.resize(len1 + len2 + len3 + 3);
        for (int i = 0; i < len1; i++) {
            if (info1.first[i].size() != 0) {
                tuple_info.first[i].resize(2);
                for (int j = 0; j < 2; j++) {
                    if (info1.first[i][j] >= len1) {
                        tuple_info.first[i][j] = info1.first[i][j] + len2 + 2;
                    } else {
                        tuple_info.first[i][j] = info1.first[i][j];
                    }
                }
            }
        }
        tuple_info.first[len1].push_back(len1 + len2 + 1);
        tuple_info.first[len1].push_back(len1 + len2 + len3 + 2);
        for (int i = 0; i < len2; i++) {
            if (info2.first[i].size() != 0) {
                tuple_info.first[len1 + 1 + i].resize(2);
                for (int j = 0; j < 2; j++) {
                    tuple_info.first[len1 + 1 + i][j] = info2.first[i][j] + len1 + 1;
                }
            }
        }
        tuple_info.first[len1 + len2 + 1].resize(2);
        tuple_info.first[len1 + len2 + 1][0] = len1;
        tuple_info.first[len1 + len2 + 1][1] = len1 + len2 + len3 + 2;
        for (int i = 0; i < len3; i++) {
            if (info1.first[i + len1].size() != 0) {
                tuple_info.first[len1 + len2 + 2 + i].resize(2);
                for (int j = 0; j < 2; j++) {
                    if (info1.first[i + len1][j] >= len1) {
                        tuple_info.first[len1 + len2 + 2 + i][j] = info1.first[len1 + i][j] + len2 + 2;
                    } else {
                        tuple_info.first[len1 + len2 + 2 + i][j] = info1.first[len1 + i][j];
                    }
                }
            }
        }
        tuple_info.first[len1 + len2 + len3 + 2].resize(2);
        tuple_info.first[len1 + len2 + len3 + 2][0] = len1;
        tuple_info.first[len1 + len2 + len3 + 2][1] = len1 + len2 + 1;

        tuple_info.second = tuple_energy(_seq[a], _seq[b], _seq[c]).first + info1.second + info2.second;

        if (len1 != 0) {
            if (!info1.first[len1 - 1].empty()) {
                tuple_info.second += stack_energy(_seq[a - 1], _seq[a]);
            }
        }
        if (len2 != 0) {
            if (!info2.first[0].empty()) {
                tuple_info.second += stack_energy(_seq[a], _seq[a + 1]);
            }
            if (!info2.first[len2 - 1].empty()) {
                tuple_info.second += stack_energy(_seq[b - 1], _seq[b]);
            }
        }
        if (len3 != 0) {
            if (!info1.first[len1].empty()) {
                tuple_info.second += stack_energy(_seq[b], _seq[b + 1]);
            }
            if (!info1.first[len1 + len3 - 1].empty()) {
                tuple_info.second += stack_energy(_seq[c - 1], _seq[c]);
            }
        }
        return tuple_info;
    }

    double stack_energy(char m, char n) {
        int a = _convert[m];
        int b = _convert[n];
        return _stack_energy(a, b);
    }

    std::pair<double, bool> pair_energy(char m, char n) {
        int a = _convert[m];
        int b = _convert[n];
        if (a == 0 && b == 1) {
            return {-484.06, true};
        } else if (a == 1 && b == 0) {
            return {-491.53, true};
        } else if (a == 2 && b == 3) {
            return {-679.05, true};
        } else if (a == 3 && b == 2) {
            return {-663.93, true};
        } else if (a == 2 && b == 1) {
            return {-200, true};
        } else if (a == 1 && b == 2) {
            return {-200, true};
        } else {
            return {-45.94, false};
        }
    }


    std::pair<double, bool> tuple_energy(char a, char b, char c) {
        auto first = pair_energy(a, b);
        auto second = pair_energy(b, c);
        return {first.first + second.first, first.second && second.second};
    }

    InfoList sub_info(const std::vector<int> &seq, double cutoff) {
        InfoList info_list;

        if (seq.size() < 3) {
            info_list.push_back(strand_info(seq));
            return info_list;
        }

        /// get best information
        auto min_en = best_info(seq)[0].second;

        /// get sub-optimal information of the first n - 1 nucleotides
        std::vector<int> temp_seq;
        std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_seq));
        auto temp_best_info = best_info(temp_seq);
        std::vector<double> ens;
        for (auto &&info: temp_best_info) {
            ens.push_back(score(info, seq.back()).second);
        }
        auto temp_en = *(std::min_element(ens.begin(), ens.end()));
        if (temp_en - min_en <= cutoff) {
            for (auto &&info: sub_info(temp_seq, cutoff - temp_en + min_en)) {
                auto temp_info = score(info, seq.back());
                if (temp_info.second - min_en <= cutoff) {
                    info_list.push_back(temp_info);
                }
            }
        }

        /// get sub-optimal information of each possible base triple
        for (int i = 0; i < seq.size() - 2; i++) {
            for (int j = i + 1; j < seq.size() - 1; j++) {
                if (seq[j] - seq[i] < _min_hairpin_size + 1 || seq.back() - seq[j] < _min_hairpin_size + 1 ||
                        tuple_energy(_seq[seq[i]], _seq[seq[j]], _seq[seq.back()]).second == false) {
                    continue;
                }

                // construct fragment 1 and fragment 2
                std::vector<int> temp_seq1, temp_seq2;
                std::copy(seq.begin(), std::next(seq.begin(), i), std::back_inserter(temp_seq1));
                std::copy(std::next(seq.begin(), j + 1), std::prev(seq.end(), 1), std::back_inserter(temp_seq1));
                std::copy(std::next(seq.begin(), i + 1), std::next(seq.begin(), j), std::back_inserter(temp_seq2));

                /// get best energy of peculiar 2D structure
                auto best_info_list1 = best_info(temp_seq1);
                auto best_info_list2 = best_info(temp_seq2);
                std::vector<double> ens;
                for (auto &&info1: best_info_list1) {
                    for (auto &&info2: best_info_list2) {
                        ens.push_back(score(info1, info2, seq[i], seq[j], seq.back()).second);
                    }
                }
                double temp_en = *(std::min_element(ens.begin(), ens.end()));

                /// push all the possible 2D structure into information list
                if (temp_en - min_en <= cutoff) {
                    auto sub_info_list1 = sub_info(temp_seq1, cutoff - temp_en + min_en);
                    auto sub_info_list2 = sub_info(temp_seq2, cutoff - temp_en + min_en);
                    for (auto &&info1: sub_info_list1) {
                        for (auto &&info2: sub_info_list2) {
                            auto new_info = score(info1, info2, seq[i], seq[j], seq.back());
                            if (new_info.second - min_en <= cutoff) {
                                info_list.push_back(new_info);
                            }
                        }
                    }
                }
            }
        }

        return info_list;
    }

};

} // namespace nuc2d
} // namespace jian

