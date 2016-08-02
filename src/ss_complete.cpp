#include <array>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <jian/nuc2d.hpp>
#include "nsp.hpp"

namespace jian {

namespace ss_complete_detail {

using pair_t = std::array<int, 2>;
using pairs_t = std::list<pair_t>;
using seq_t = std::string;
using ss_t = std::string;

void set_pairs(pairs_t & pairs, const ss_t & ss) {
    auto & keys = NucSS::instance().paired_keys;
    std::vector<std::deque<int>> v(keys.size());
    int i = 0;
    for (auto && c : ss) {
        auto it1 = std::find_if(keys.begin(), keys.end(), [&c](auto && key){return key.first == c;});
        if (it1 != keys.end()) {
            int n = std::distance(keys.begin(), it1);
            v[n].push_back(i);
        } else {
            auto it2 = std::find_if(keys.begin(), keys.end(), [&c](auto && key){return key.second == c;});
            if (it2 != keys.end()) {
                int n = std::distance(keys.begin(), it2);
                pairs.push_back({v[n].back(), i});
                v[n].pop_back();
            } else {
            }
        }
        i++;
    }
}

void print_pairs(const pairs_t & pairs) {
    for (auto && pair : pairs) {
        std::cout << pair[0]+1 << ' ' << pair[1]+1 << std::endl;
    }
}

bool is_paired(const char & a, const char & b) {
    return (a == 'A' && b == 'U') ||
           (a == 'U' && b == 'A') ||
           (a == 'G' && b == 'C') ||
           (a == 'C' && b == 'G') ||
           (a == 'G' && b == 'U') ||
           (a == 'U' && b == 'G');
}

void complete(pairs_t & pairs, const seq_t & seq) {
    int l = seq.size();
    auto it = pairs.begin();
    while (true) {
        if (it != pairs.end()) {
            auto & pair = *it;
            if (!is_paired(seq[pair[0]], seq[pair[1]])) {
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

}

REGISTER_NSP_COMPONENT(ss_complete) {
    using namespace ss_complete_detail;
    seq_t seq = par["seq"][0];
    ss_t ss = par["ss"][0];
    pairs_t pairs;
    set_pairs(pairs, ss);
    complete(pairs, seq);
    print_pairs(pairs);
}

} // namespace jian

