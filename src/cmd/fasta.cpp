#include "nsp.hpp"
#include "fasta.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(fasta) {
    auto g = par.getv("global");

    Str fasta_file = g[1];
    Fasta fa;
    fasta_read(fa, fasta_file);

    if (par.has("random", "rand")) {
        std::random_shuffle(fa.begin(), fa.end());
    }

    Int n = JN_INT(par.get("n", "num"));
    fasta_write(JN_OUT, fa, n);
}

}

