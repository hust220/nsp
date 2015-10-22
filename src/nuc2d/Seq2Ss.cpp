#include "Seq2Ss.h"

namespace jian {
namespace nuc2d {

std::vector<std::pair<int, int>> Seq2Ss::operator ()(std::string seq) {
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
    backtrack(0, _len - 1);
    for (auto &&pair: _pairs) {
        std::cout << pair.first << ' ' << pair.second << std::endl;
    }
    return _pairs;
}

int Seq2Ss::mms(int m, int n) {
    if (n <= m + _min_hairpin_size) {
        return 0;
    }
    int max_size = get_mms(m, n - 1);
    for (int i = m; i < n - _min_hairpin_size; i++) {
        int temp = get_mms(m, i - 1) + get_mms(i + 1, n - 1) + score(i, n);
        if (temp > max_size) {
            max_size = temp;
        }
    }
    return max_size;
}

int Seq2Ss::score(int m, int n) {
    int a = _types[m] + _types[n];
    int b = _types[m] * _types[n];
    if (a == 1 || a == 5 || b == 2) {
        return 1;
    } else {
        return 0;
    }
}

void Seq2Ss::backtrack(int m, int n) {
    if (n <= m + _min_hairpin_size) {
        return;
    }
    if (get_mms(m, n - 1) == get_mms(m, n)) {
        backtrack(m, n - 1);
    } else {
        for (int i = m; i < n - _min_hairpin_size; i++) {
            int size = get_mms(m, i - 1) + get_mms(i + 1, n - 1) + score(i, n);
            if (size == get_mms(m, n)) {
                _pairs.push_back(make_pair(i, n));
                backtrack(m, i - 1);
                backtrack(i + 1, n - 1);
                break;
            }
        }
    }
}

int Seq2Ss::get_mms(int m, int n) {
    if (n <= m + _min_hairpin_size) {
        return 0;
    } else {
        return _mms(m, n);
    }
    
}
        
} /// namespace nuc2d
} /// namespace jian
