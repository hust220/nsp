#include "nsp.hpp"
#include "pdb.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(seq_ident) {
    auto g = par.getv("global");
    Str seq1 = g[1];
    Str seq2 = g[2];
    Int l1 = size(seq1);
    Int l2 = size(seq2);
    if (l1 != l2) throw "Length of two sequences should be equal!";
    Int n = 0;
    for (Int i = 0; i < l1; i++) {
        if (seq1[i] == seq2[i]) n++;
    }
    JN_OUT << Num(n) / l1 << std::endl;
}

}

