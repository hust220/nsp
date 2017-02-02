#include <nsp/lrsp.hpp>
#include <nsp/dca.hpp>
#include <jian/utils/log.hpp>
#include "nsp.hpp"
#include <jian/utils/traits.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(ssp_dca) {
    S seq = par["seq"][0];
    double k = std::stof(par["k"][0]);
    S di_file = par["di"][0];

    LOGI << "Sequence:" << std::endl;
    LOGI << seq << std::endl;

    int size = seq.size();
    if (par.has("k")) {
        size = int(size * k);
    }
    LOGI << "Read first " << size << " pairs:" << std::endl;
    auto && pairs = dca::pairs_from_file(di_file, size);

    LOGI << "Print pairs:" << std::endl;
    dca::print_pairs(JN_OUT, pairs);

    LOGI << "Converting DI to ss: " << std::endl;
    dca::ss_t ss = dca::pairs_to_ss(pairs, seq.size());
    pairs = dca::pairs_from_ss(ss);
    LOGI << ss << std::endl;

    LOGI << "Extending base pairs: " << std::endl;
    lrsp::ss_complete(pairs, seq);
    ss = dca::pairs_to_ss(pairs, seq.size());
    LOGI << ss << std::endl;

    //LOGI << "Removing small hairpins: " << std::endl;
    //ss = lrsp::ss_del_sm_hp(ss);
    //LOGI << ss << std::endl;

    LOGI << "SS Prediction: " << std::endl;
    ss = lrsp::ss_hp_pred(seq, ss);
    LOGI << ss << std::endl;
}

END_JN

