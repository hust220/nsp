#include "nsp.hpp"
#include "fasta.hpp"

BEGIN_JN

static void msa_trim(Fasta &fa) {
    Int a = -1, b = -1;
    for (auto && it : fa) {
        if (it.name == "TARGET") {
            for (Int i = 1; i < size(it.seq); i++) {
                if (a == -1 && it.seq[i-1] == '-' && it.seq[i] != '-') {
                    a = i;
                }
                if (it.seq[i] != '-') {
                    b = i;
                }
            }
        }
    }
    for (auto && it : fa) {
        it.seq = it.seq.substr(a, b-a+1);
    }
}

REGISTER_NSP_COMPONENT(msa) {
    auto g = par.getv("global");

    Str cmd = g[1];

    if (cmd == "trim") {
        Str fasta_file = g[2];
        Fasta fa;
        fasta_read(fa, fasta_file);
        msa_trim(fa);
        fasta_write(JN_OUT, fa);
    }

}

END_JN

