#include <jian/lrsp.hpp>
#include <jian/dca.hpp>
#include <jian/utils/log.hpp>
#include "nsp.hpp"

namespace jian {

REGISTER_NSP_COMPONENT(ss_dca) {
    LOGI << "Sequence:" << std::endl;
    dca::seq_t seq = par["seq"][0];
    LOGI << seq << std::endl;

    int size = seq.size();
    if (par.has("k")) {
        size *= std::stof(par["k"][0]);
    }
    auto && pairs = dca::pairs_from_file(par["di"][0], size);

    LOGI << "Converting DI to ss: " << std::endl;
    dca::ss_t ss = dca::pairs_to_ss(pairs, seq.size());
    pairs = dca::pairs_from_ss(ss);
    LOGI << ss << std::endl;

    LOGI << "Extending base pairs: " << std::endl;
    lrsp::ss_complete(pairs, seq);
    ss = dca::pairs_to_ss(pairs, seq.size());
    LOGI << ss << std::endl;

    LOGI << "Removing small hairpins: " << std::endl;
    ss = lrsp::ss_del_sm_hp(ss);
    LOGI << ss << std::endl;

    LOGI << "SS Prediction: " << std::endl;
    ss = lrsp::ss_hp_pred(seq, ss);
    LOGI << ss << std::endl;
}

} // namespace jian

