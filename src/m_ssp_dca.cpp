#include <nsp/lrsp.hpp>
#include <nsp/dca.hpp>
#include <jian/utils/log.hpp>
#include "nsp.hpp"
#include <jian/utils/traits.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(ssp_dca) {
    Num k = 0.3;
    par.set(k, "k");

    Str seq = par.get("seq");
    Str di_file = par.get("di");

    LOGI << "Sequence:" << std::endl;
    LOGI << seq << std::endl;

    int size = seq.size() * k;
    LOGI << "Read first " << size << " pairs:" << std::endl;
    auto && pairs = dca::pairs_from_file(di_file, size);

    LOGI << "Print pairs:" << std::endl;
    dca::print_pairs(std::cout, pairs);

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

    JN_OUT << ss << std::endl;
}

END_JN

