#include "nsp.hpp"
#include <jian/utils/file.hpp>
#include <jian/matrix.hpp>
#include <jian/scoring/score_psb.hpp>
#include <jian/cg.hpp>
#include <jian/pdb.hpp>

namespace jian {

static void score(const std::string &name) {
    Chain chain = CGpsb::chain(read_model_to_chain(name));
    int len = chain.size();
    double e = 0, s;
    for (int i = 0; i < len - 1; i++) {
        s = scoring::score_stacking_psb(chain[i], chain[i+1]);
        e += s;
        for (int j = i + 2; j < len; j++) {
            s = scoring::score_pairing_psb(chain[i], chain[j]);
            e += s;
        }
    }
    std::cout << e << std::endl;
}

static void score_res(const std::string &name) {
    Chain chain = CGpsb::chain(read_model_to_chain(name));
    int len = chain.size();
    Mat m = Mat::Zero(len, len);
    double s;
    for (int i = 0; i < len - 1; i++) {
        s = scoring::score_stacking_psb(chain[i], chain[i+1]);
        m(i, i+1) = s;
        for (int j = i + 2; j < len; j++) {
            s = scoring::score_pairing_psb(chain[i], chain[j]);
            m(i, j) = s;
        }
    }
    std::cout << m << std::endl;
}

REGISTER_NSP_COMPONENT(score_psb) {
    if (par.has("s")) {
        std::string name = par.get("s");
        if (par.has("res")) {
            score_res(name);
        } else {
            score(name);
        }
    } else if (par.has("list")) {
        EACH_SPLIT_LINE(par["list"][0].c_str(), " ",
            score(F[0]);
        );
    }
}

REGISTER_NSP_COMPONENT(train_psb) {
    scoring::train_psb(par["list"][0]);
}

REGISTER_NSP_COMPONENT(freqs_psb) {
    scoring::print_freqs_psb();
}

} // namespace jian
















