#include "../dca.hpp"
#include <array>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include "../nuc2d.hpp"
#include "ss_complete.hpp"

BEGIN_JN
namespace lrsp {

std::map<char, int> index_base {{'X', 0}, {'A', 1}, {'U', 2}, {'T', 2}, {'G', 4}, {'C', 8}};

bool is_paired(const char & a, const char & b) {
    int n = index_base[a] | index_base[b];
    return n == 3 || n == 12 || n == 6;
}

void ss_complete(pairs_t & pairs, const seq_t & seq) {
    int l = seq.size();
    auto it = pairs.begin();
    while (true) {
        if (it != pairs.end()) {
            auto & pair = *it;
            if (pair[1] - pair[0] - 1 <= 3 || !is_paired(seq[pair[0]], seq[pair[1]])) {
                it = pairs.erase(it);
                continue;
            }
            pair_t p1 {pair[0] - 1, pair[1] + 1};
            pair_t p2 {pair[0] + 1, pair[1] - 1};
            int flag = 0;
            if (p1[0] >= 0 && p1[1] < l && is_paired(seq[p1[0]], seq[p1[1]])) {
                if (std::find(pairs.begin(), pairs.end(), p1) == pairs.end()) {
                    pairs.push_back(p1);
                }
                flag++;
            }
            if (p2[1] - p2[0] > 3 && is_paired(seq[p2[0]], seq[p2[1]])) {
                if (std::find(pairs.begin(), pairs.end(), p2) == pairs.end()) {
                    pairs.push_back(p2);
                }
                flag++;
            }
            if (flag == 0) {
                it = pairs.erase(it);
            } else {
                it++;
            }
        } else {
            break;
        }
    }
}

S ss_complete(S seq, S ss) {
    pairs_t pairs = dca::pairs_from_ss(ss);
    ss_complete(pairs, seq);
    return dca::pairs_to_ss(pairs, seq.size());
}

} // namespace lrsp
END_JN

