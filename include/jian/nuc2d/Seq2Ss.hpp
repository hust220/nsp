#ifndef JIAN_NUC2D_SEQ2SS_H
#define JIAN_NUC2D_SEQ2SS_H

#include "../etl.hpp"

namespace jian {
namespace nuc2d {

class Seq2Ss {
public:
    using Pair = std::pair<int, int>;
    using Pairs = std::vector<Pair>;
    using PairList = std::vector<int>;
    using PairInfo = std::pair<PairList, double>;
    using InfoList = std::vector<PairInfo>;

    std::map<int, InfoList> _info;

    MatrixXf _mms;
    std::string _seq;
    std::vector<int> _types;
    int _len;
    int _min_hairpin_size;
    double _cutoff;
    Matrix4f _pair_energy;
    Matrix4f _stack_energy;

    Seq2Ss() {
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
        std::map<char, int> temp_map{{'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3}};
        std::transform(_seq.begin(), _seq.end(), std::back_inserter(_types), [&](char c){return temp_map[c];});

        auto info = best_info(0, len - 1);
        auto info_list = sub_info(0, len - 1, -_cutoff * info[0].second);
        std::sort(info_list.begin(), info_list.end(), [](const PairInfo &info1, const PairInfo &info2){return info1.second < info2.second;});
        for (auto &&info: info_list) {
            for (int i = 0; i < info.first.size(); i++) {
                if (info.first[i] == -1) {
                    std::cout << '.';
                } else if (info.first[i] > i) {
                    std::cout << '(';
                } else {
                    std::cout << ')';
                }
                std::cout << ' ';
            }
            std::cout << info.second << std::endl;
        }
        return info_list;
    }

    InfoList best_info(int m, int n) {
        InfoList info_list;
        int len = _seq.size();

        if (n <= m + _min_hairpin_size) {
            info_list.push_back(strand_info(m, n));
            return info_list;
        }
        if (_info.count(m * len + n)) {
            return _info[m * len + n];
        }
        InfoList temp_info_list;
        std::vector<int> temp_vec;
        auto best_info_list = best_info(m, n - 1);
        for (auto &&info: best_info_list) {
            temp_info_list.push_back(score(info, n));
        }
        for (int i = m; i < n - _min_hairpin_size; i++) {
            if (pair_energy(i, n).second == false) continue;
            auto info_list1 = best_info(m, i - 1);
            auto info_list2 = best_info(i + 1, n - 1);
            for (auto &&info1: info_list1) {
                for (auto &&info2: info_list2) {
                    temp_info_list.push_back(score(info1, i, info2, n));
                }
            }
        }
        temp_vec.push_back(0);
        for (int i = 1; i < temp_info_list.size(); i++) {
            if (temp_info_list[i].first.size() == 0) continue;
            double en1 = temp_info_list[i].second;
            double en2 = temp_info_list[temp_vec[0]].second;
            if (en1 < en2) {
                temp_vec.clear();
                temp_vec.push_back(i);
            } else if (en1 == en2) {
                temp_vec.push_back(i);
            }
        }
        for (auto &&i: temp_vec) {
            info_list.push_back(temp_info_list[i]);
        }
        _info[m * len + n] = info_list;
        return info_list;
    }

    PairInfo strand_info(int m, int n) {
        PairInfo pair_info;
        if (m > n) {
            return pair_info;
        } else {
            int len = n - m + 1;
            pair_info.first.resize(len);
            for (int i = 0; i < len; i++) {
                pair_info.first[i] = -1;
            }
            pair_info.second = 0;
            return pair_info;
        }
    }

    PairInfo score(const PairInfo &info, int m) {
        int len = info.first.size() + 1;
        PairInfo pair_info;
        pair_info.first.resize(len);
        std::copy(info.first.begin(), info.first.end(), pair_info.first.begin());
        pair_info.first[len - 1] = -1;
        pair_info.second = info.second + 0;
        return pair_info;
    }

    PairInfo score(const PairInfo &info1, int m, const PairInfo &info2, int n) {
        int len1 = info1.first.size();
        int len2 = info2.first.size();

        PairInfo pair_info;
        pair_info.first.resize(len1 + len2 + 2);
        for (int i = 0; i < len1; i++) pair_info.first[i] = info1.first[i];
        for (int i = 0; i < len2; i++) pair_info.first[i + len1 + 1] = (info2.first[i] == -1 ? -1 : info2.first[i] + len1 + 1);
        pair_info.first[len1] = len1 + 1 + len2;
        pair_info.first[len1 + 1 + len2] = len1;
        pair_info.second = pair_energy(m, n).first;
        if (info1.first.size() != 0) {
            pair_info.second += info1.second;
            if (info1.first.back() != -1) {
                pair_info.second += stack_energy(m - 1, m);
            }
        }
        if (info2.first.size() != 0) {
            pair_info.second += info2.second;
            if (info2.first.front() != -1) {
                pair_info.second += stack_energy(m, m + 1);
            }
            if (info2.first.back() != -1) {
                pair_info.second += stack_energy(n - 1, n);
            }
        }
        return pair_info;
    }

    double stack_energy(int m, int n) {
        int a = _types[m];
        int b = _types[n];
        return _stack_energy(a, b);
    }

    std::pair<double, bool> pair_energy(int m, int n) {
        int a = _types[m];
        int b = _types[n];
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

    InfoList sub_info(int m, int n, double cutoff) {
        InfoList temp_vec;
        if (n <= m + _min_hairpin_size) {
            temp_vec.push_back(strand_info(m, n));
            return temp_vec;
        }

        auto min_en = best_info(m, n)[0].second;
        auto temp_info_list = best_info(m, n - 1);
        std::vector<double> ens;
        for (auto &&info: temp_info_list) {
            ens.push_back(score(info, n).second);
        }
        auto temp_en = *(std::min_element(ens.begin(), ens.end()));
        if (temp_en - min_en <= cutoff) {
            auto temp_info_list = sub_info(m, n - 1, cutoff - temp_en + min_en);
            for (auto &&info: temp_info_list) {
                temp_vec.push_back(score(info, n));
            }
        }
        for (int i = m; i < n - _min_hairpin_size; i++) {
            if (pair_energy(i, n).second == false) continue;
            auto best_info_list1 = best_info(m, i - 1);
            auto best_info_list2 = best_info(i + 1, n - 1);
            std::vector<double> ens;
            for (auto &&info1: best_info_list1) {
                for (auto &&info2: best_info_list2) {
                    ens.push_back(score(info1, i, info2, n).second);
                }
            }
            auto temp_en = *(std::min_element(ens.begin(), ens.end()));
            if (temp_en - min_en <= cutoff) {
                auto info_list1 = sub_info(m, i - 1, cutoff - temp_en + min_en);
                auto info_list2 = sub_info(i + 1, n - 1, cutoff - temp_en + min_en);
                for (int j = 0; j < info_list1.size(); j++) {
                    for (int k = 0; k < info_list2.size(); k++) {
                        auto temp_info2 = score(info_list1[j], i, info_list2[k], n);
                        if (temp_info2.second - min_en <= cutoff) {
                            temp_vec.push_back(temp_info2);
                        }
                    }
                }
            }
        }
        return temp_vec;
    }

};

} /// namespace nuc2d
} /// namespace jian

#endif
