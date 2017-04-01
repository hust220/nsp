#include "bires.hpp"
#include "bires_models.hpp"
#include "../nuc2d/NASS.hpp"
#include <jian/utils/rand.hpp>

BEGIN_JN

static Int seq_code(Char seq1, Char seq2) {
    // fill
    static Map<Char, Int> m{{'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3}};
    return m[seq1]*4+m[seq2];
}

static Pair<Int, Int> ss_code(Char c) {
    const auto & pks = NASS::instance().paired_keys;
    const auto & uks = NASS::instance().unpaired_keys;
    auto it = std::find_if(pks.begin(), pks.end(), [&c](auto && k){return k.first == c;});
    if (it != pks.end()) return {1, std::distance(pks.begin(), it)};
    else {
        auto it2 = std::find_if(pks.begin(), pks.end(), [&c](auto && k){return k.second == c;});
        if (it2 != pks.end()) return {2, std::distance(pks.begin(), it2)};
        else return {0, 0};
    }
}

static Int ss_code(Char ss1, Char ss2) {
    auto c1 = ss_code(ss1);
    auto c2 = ss_code(ss2);
    if      (c1.first == 1 && c2.first == 1 && c1.second == c2.second) return 0;
    else if (c1.first == 1 && c2.first == 2 && c1.second == c2.second) return 1;
    else if (c1.first == 2 && c2.first == 2 && c1.second == c2.second) return 2;
    else return 3;
}

Int bires_code(Char seq1, Char seq2, Char ss1, Char ss2) {
    return seq_code(seq1, seq2) * 4 + ss_code(ss1, ss2);
}

const Chain &bires_chain(Char seq1, Char seq2, Char ss1, Char ss2) {
    auto & ms = bires_models();
    auto & m = ms[bires_code(seq1, seq2, ss1, ss2)];
    Int l = size(m);
    /*
    if (l == 0) {
        std::cout << bires_code(seq1, seq2, ss1, ss2) << std::endl;
    }
    */
    Int n = Int(JN_ rand() * l);
    return m[n];
}


END_JN


