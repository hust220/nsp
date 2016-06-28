#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>
#include <string>
#include <unordered_map>
#include "../utils/ls.hpp"
#include "../matrix.hpp"

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

class PredTriSS {
public:
    using tuple_t = std::vector<int>;
    using tuples_t = std::vector<tuple_t>;
    using info_t = std::pair<tuples_t, double>;
    using infos_t = std::vector<info_t>;
    using seq_t = std::vector<int>;

    std::unordered_map<std::vector<int>, infos_t, seq2tri::MyHash> _info;

    std::string _seq;
    std::vector<int> _types;
    std::vector<std::tuple<int, int, int>> _pairs;
    std::map<char, int> _convert{{'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3}};
    std::hash<std::string> hash_fn;
    int _min_hairpin_size;
    double _cutoff;
    Eigen::Matrix4f _pair_energy;
    Eigen::Matrix4f _stack_energy;

    PredTriSS() {
        _cutoff = 0.4;
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

    void run(std::string seq) {
        _seq = seq;
        seq_t vec(_seq.size());
        std::iota(vec.begin(), vec.end(), 0);
        auto &&best_infos = best_info(vec);
        auto &&sub_infos = sub_info(vec, -_cutoff * best_infos[0].second);
        std::deque<std::string> ss;
        for (auto && i : sub_infos) ss.push_back(dbn(i.first, _seq.size()));
        std::deque<double> sc;
        for (auto && i : ss) sc.push_back(score_ss(i));
        std::vector<int> arr(sub_infos.size());
        std::iota(arr.begin(), arr.end(), 0);
        std::sort(arr.begin(), arr.end(), [&](int i, int j){return sc[i] < sc[j];});
        for (auto it = arr.begin(); it < arr.begin() + 10 && it < arr.end(); it++) {
            std::cout << ss[*it] << ' ' << sc[*it] << std::endl;
        }
    }

    double score_ss(const std::string &ss) {
        seq_t loop;
        double e = 0;
        for (int i = 0; i < ss.size(); i++) {
            if (i < ss.size() - 1) {
                if (ss[i] != '0' && ss[i+1] != '0') {
                    e += -1;
                }
            }
            if (ss[i] == '0') {
                loop.push_back(i);
            } else {
                if (ss[i] == '1') e += -1;
                if (!loop.empty()) {
                    e += 3 + std::log(loop.size());
                    loop.clear();
                }
            }
        }
        if (!loop.empty()) {
            e += 1 + std::log(loop.size());
            loop.clear();
        }
        return e;
    }

    std::string dbn(const tuples_t &tuples, int len) {
        auto init_ss = [&](auto &&ss){
            for (int i = 0; i < tuples.size(); i++) {
                if (tuples[i].size() != 0) {
                    std::vector<int> v {tuples[i][0], tuples[i][1], tuples[i][2]};
                    std::sort(v.begin(), v.end());
                    ss[v[0]] = '1'; ss[v[1]] = '2'; ss[v[2]] = '3';
                }
            }
        };

        auto count_purine = [&](auto &&ss, auto &&v){
            for (int i = 0; i < ss.size(); i++) {
                if (_seq[i] == 'A' || _seq[i] == 'G') {
                    if (ss[i] == '1') v[0]++;
                    else if (ss[i] == '2') v[1]++;
                    else if (ss[i] == '3') v[2]++;
                }
            }
        };

        auto reset_ss = [](auto &&ss, auto &&m) {
            for (auto && c : ss) {
                if (c == '1') c = m[0];
                else if (c == '2') c = m[1];
                else if (c == '3') c = m[2];
            }
        };

        std::string ss(len, '0');
        init_ss(ss);
        std::array<int, 3> v {0, 0, 0};
        count_purine(ss, v);
        std::array<int, 3> arr {0, 1, 2};
        std::sort(arr.begin(), arr.end(), [&v](int a, int b){return v[a] > v[b];});
        if (arr[1] == 1) std::swap(arr[0], arr[1]);
        std::array<char, 3> m;
        m[arr[0]] = '1'; m[arr[1]] = '3'; m[arr[2]] = '2';
        reset_ss(ss, m);
        return ss;
    }

    infos_t best_info(const seq_t &seq) {
        infos_t info_list;

        if (seq.size() < 3) {
            info_list.push_back(strand_info(seq));
            return info_list;
        }

        if (_info.count(seq)) {
            return _info[seq];
        }

        /// add all the best information of first n - 1 nucleotides
        infos_t temp_info_list;
        std::vector<int> temp_vec;
        std::copy(seq.begin(), std::prev(seq.end(), 1), std::back_inserter(temp_vec));
        auto temp_best_info = best_info(temp_vec);
        for (auto &&info: temp_best_info) {
            temp_info_list.push_back(score(seq, info));
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

                // calculate best information
                auto best_info1 = best_info(temp_seq1);
                auto best_info2 = best_info(temp_seq2);
                for (auto &&info1: best_info1) {
                    for (auto &&info2: best_info2) {
                        temp_info_list.push_back(score(seq, info1, info2, seq[i], seq[j], seq.back()));
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

    info_t strand_info(const seq_t & seq) {
        info_t tuple_info;
        if (seq.empty()) {
            tuple_info.second = 0;
        } else {
            tuple_info.second = std::log(seq.size()) + 1;
        }
        return tuple_info;
    }

    info_t score(const seq_t & seq, info_t info) {
        return info;
    }

    info_t score(const seq_t & seq, const info_t &info1, const info_t &info2, int a, int b, int c) {
        info_t info;
        for (auto && i : info1.first) info.first.push_back(i);
        for (auto && i : info2.first) info.first.push_back(i);
        info.second = info1.second + info2.second;
        if (info.first.empty() || ({auto &t = info.first.back(); a > t[0] && b < t[1] && c > t[2];})) {
            info.second += -1;
            info.first.push_back({a, b, c});
        }
        return info;
    }

    double stack_energy(char m, char n) {
        int a = _convert[m];
        int b = _convert[n];
        return _stack_energy(a, b);
    }

    std::pair<double, bool> tuple_energy(char a, char b, char c) {
        int x = _convert[a], y = _convert[b], z = _convert[c];
        auto att = [](int n)->std::pair<double,bool>{ if (n == 1) return {-484.06*1.5, true}; else return {0, false}; };
        auto gcc = [](int n)->std::pair<double,bool>{ if (n == 2 || n == 3) return {-679.05*1.5, true}; else return {0, false}; };
        if (x + y == 1) return att(z);
        else if (x + y == 5) return gcc(z);
        else if (y + z == 1) return att(x);
        else if (y + z == 5) return gcc(x);
        else return {0, false};
    }

    infos_t sub_info(const seq_t &seq, double cutoff) {
        infos_t info_list;

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
            ens.push_back(score(seq, info).second);
        }
        auto temp_en = *(std::min_element(ens.begin(), ens.end()));
        if (temp_en - min_en <= cutoff) {
            for (auto &&info: sub_info(temp_seq, cutoff - temp_en + min_en)) {
                auto temp_info = score(seq, info);
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
                        ens.push_back(score(seq, info1, info2, seq[i], seq[j], seq.back()).second);
                    }
                }
                double temp_en = *(std::min_element(ens.begin(), ens.end()));

                /// push all the possible 2D structure into information list
                if (temp_en - min_en <= cutoff) {
                    auto sub_info_list1 = sub_info(temp_seq1, cutoff - temp_en + min_en);
                    auto sub_info_list2 = sub_info(temp_seq2, cutoff - temp_en + min_en);
                    for (auto &&info1: sub_info_list1) {
                        for (auto &&info2: sub_info_list2) {
                            auto new_info = score(seq, info1, info2, seq[i], seq[j], seq.back());
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

