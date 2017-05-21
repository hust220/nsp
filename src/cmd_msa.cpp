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

static void msa_del_gap(Fasta &fa) {
    Map<Int, Bool> m;
    for (auto && it : fa) {
        if (it.name == "TARGET") {
            for (Int i = 0; i < size(it.seq); i++) {
                if (it.seq[i] == 'A' || it.seq[i] == 'U' || it.seq[i] == 'G' || it.seq[i] == 'C') {
                    m[i] = true;
                }
                else {
                    m[i] = false;
                }
            }
        }
    }
    for (auto && it : fa) {
        std::stringstream stream;
        for (Int i = 0; i < size(it.seq); i++) {
            if (m[i]) {
                stream << it.seq[i];
            }
        }
        it.seq = stream.str();
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
        fasta_write(JN_OUT, fa, -1, -1);
    }
    else if (cmd == "del_gap") {
        Str fasta_file = g[2];
        Fasta fa;
        fasta_read(fa, fasta_file);
        msa_del_gap(fa);
        fasta_write(JN_OUT, fa, -1, -1);
    }

}

END_JN

