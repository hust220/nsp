#include "nsp.hpp"
#include "pdb.hpp"
#include "g4dna_ssp.hpp"

BEGIN_JN

static void handle_g4dna_ssp(Str seq) {
    auto infos = g4dna_ssp(seq);
    for (auto && info : infos) {
        for (auto && frag : info) {
            Int i = 0; for (; i+1 < size(frag); i++) JN_OUT << frag[i] << '-'; JN_OUT << frag[i] << ' ';
        }
        JN_OUT << std::endl;
    }
}

REGISTER_NSP_COMPONENT(g4dna) {
    auto g = par.getv("global");
    Str t = g[1];
    if (t == "ssp") {
        Str seq = g[2];
        handle_g4dna_ssp(seq);
    }
}

END_JN

