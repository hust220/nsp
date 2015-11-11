#include "Seq2Ss.h"

namespace jian {
namespace nuc2d {

Seq2Ss::Seq2Ss() {
	_pair_energy << -45.94, -484.06, -45.94, -45.94,
			        -491.53, -45.94, -200, -45.94,
					-45.94, -200, -45.94, -679.05,
					-45.94, -45.94, -663.93, -45.94;
	_stack_energy << -625.57, -578.22, -772.93, -650.16,
			         -700.71, -580.50, -753.00, -525.07,
					 -742.75, -700.00, -837.22, -810.87,
					 -642.19, -735.33, -851.11, -509.75;
}

void Seq2Ss::operator ()(std::string seq) {
    _seq = seq;
    int len = seq.size();
    std::map<char, int> temp_map{{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
    std::transform(_seq.begin(), _seq.end(), std::back_inserter(_types), [&](char c){return temp_map[c];});

    auto info = best_info(0, len - 1);
    auto info_list = sub_info(0, len - 1, -_cutoff * info.second);
    std::sort(info_list.begin(), info_list.end(), [](const PairInfo &info1, const PairInfo &info2){return info1.second < info2.second;});
    for (auto &&info: info_list) {
    	for (auto i: info.first) {
    		std::cout << i << ' ';
    	}
    	std::cout << info.second << std::endl;
    }
}

Seq2Ss::PairInfo Seq2Ss::best_info(int m, int n) {
    PairInfo pair_info;
    int len = _seq.size();
    if (n <= m + _min_hairpin_size) {
        return strand_info(m, n);
    }
    if (_info.count(m * len + n)) {
        return _info[m * len + n];
    }
    auto min_score = score(best_info(m, n - 1), n);
    for (int i = m; i < n - _min_hairpin_size; i++) {
        auto temp = score(best_info(m, i - 1), i, best_info(i + 1, n - 1), n);
        if (min_score.first.size() != 0 && temp.first.size() != 0 && temp.second < min_score.second) {
            min_score = temp;
        }
    }
    _info[m * len + n] = min_score;
    return min_score;
}

Seq2Ss::PairInfo Seq2Ss::strand_info(int m, int n) {
    PairInfo pair_info;
    if (m > n) {
        return pair_info;
    } else {
        int len = n - m + 1;
        pair_info.first.resize(len);
        for (int i = 0; i < len; i++) {
            pair_info.first[i] = -1;
        }
        pair_info.second = len;
        return pair_info;
    }
}

Seq2Ss::PairInfo Seq2Ss::score(const Seq2Ss::PairInfo &info, int m) {
	int len = info.first.size() + 1;
	PairInfo pair_info;
	pair_info.first.resize(len);
	std::copy(info.first.begin(), info.first.end(), pair_info.first.begin());
	pair_info.first[len - 1] = -1;
	pair_info.second = info.second + 0;
	return pair_info;
}

Seq2Ss::PairInfo Seq2Ss::score(const Seq2Ss::PairInfo &info1, int m, const Seq2Ss::PairInfo &info2, int n) {
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

double Seq2Ss::stack_energy(int m, int n) {
	int a = _types[m];
	int b = _types[n];
	return _stack_energy(a, b);
}

std::pair<double, bool> Seq2Ss::pair_energy(int m, int n) {
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

Seq2Ss::InfoList Seq2Ss::sub_info(int m, int n, double cutoff) {
    InfoList temp_vec;
    if (n <= m + _min_hairpin_size) {
        temp_vec.push_back(strand_info(m, n));
        return temp_vec;
    }

    auto min_info = best_info(m, n);
    auto temp_info = best_info(m, n - 1);
    if (temp_info.second - min_info.second <= cutoff) {
        auto temp_info_list = sub_info(m, n - 1, cutoff - temp_info.second + min_info.second);
        for (auto &&info: temp_info_list) {
        	temp_vec.push_back(score(info, n));
        }
    }
    for (int i = m; i < n - _min_hairpin_size; i++) {
    	auto info1 = best_info(m, i - 1);
    	auto info2 = best_info(i + 1, n - 1);
    	auto temp_info = score(info1, i, info2, n);
        if (temp_info.second - min_info.second <= cutoff) {
            auto info_list1 = sub_info(m, i - 1, cutoff - temp_info.second + min_info.second);
            auto info_list2 = sub_info(i + 1, n - 1, cutoff - temp_info.second + min_info.second);
            for (int j = 0; j < info_list1.size(); j++) {
            	for (int k = 0; k < info_list2.size(); k++) {
            		auto temp_info2 = score(info_list1[j], i, info_list2[k], n);
                    if (temp_info2.second - min_info.second <= cutoff) {
                    	temp_vec.push_back(temp_info2);
                    }
            	}
            }
        }
    }
    return temp_vec;
}

} /// namespace nuc2d
} /// namespace jian
