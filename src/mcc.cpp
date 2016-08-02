#include "nsp.hpp"
#include <jian/nuc2d.hpp>

namespace jian {

namespace mcc_detail {

using ss_t = std::string;
using pair_t = std::pair<int, int>;
using pairs_t = std::vector<pair_t>;

// interaction network type
struct in_t {
    double tp, fp, fn;
};

pairs_t get_pairs(const ss_t &ss) {
    auto & v = NucSS::instance().paired_keys;
    std::vector<std::vector<int>> st(v.size());
    pairs_t pairs;
    int i = 0;
    for (auto && c : ss) {
        if (c == '&') continue;
        auto it_l = std::find_if(v.begin(), v.end(), [&c](auto && pair){return pair.first == c;});
        auto it_r = std::find_if(v.begin(), v.end(), [&c](auto && pair){return pair.second == c;});
        if (it_l != v.end()) {
            int n = std::distance(v.begin(), it_l);
            st[n].push_back(i);
        } else if (it_r != v.end()) {
            int n = std::distance(v.begin(), it_r);
            pairs.push_back({st[n].back(), i});
            st[n].pop_back();
        }
        i++;
    }
    return pairs;
}

in_t get_in(const ss_t &nat, const ss_t &pred) {
    in_t in;
    pairs_t pairs_nat = get_pairs(nat);
    pairs_t pairs_pred = get_pairs(pred);
    in.tp = std::count_if(pairs_pred.begin(), pairs_pred.end(), [&pairs_nat](auto && pair){
        return std::find(pairs_nat.begin(), pairs_nat.end(), pair) != pairs_nat.end();
    });
    in.fp = pairs_pred.size() - in.tp;
    in.fn = std::count_if(pairs_nat.begin(), pairs_nat.end(), [&pairs_pred](auto && pair){
        return std::find(pairs_pred.begin(), pairs_pred.end(), pair) == pairs_pred.end();
    });
    return in;
}

double sty(const in_t &in) {
    return in.tp / (in.tp + in.fn);
}

double ppv(const in_t &in) {
    return in.tp / (in.tp + in.fp);
}

double mcc(const in_t &in) {
    return std::sqrt(sty(in) * ppv(in));
}

} // namespace mcc_detail

REGISTER_NSP_COMPONENT(mcc) {
    mcc_detail::ss_t nat = par["nat"][0];
    mcc_detail::ss_t pred = par["pred"][0];
    std::cout << mcc_detail::mcc(mcc_detail::get_in(nat, pred)) << std::endl;
}

REGISTER_NSP_COMPONENT(sty) {
    mcc_detail::ss_t nat = par["nat"][0];
    mcc_detail::ss_t pred = par["pred"][0];
    std::cout << mcc_detail::sty(mcc_detail::get_in(nat, pred)) << std::endl;
}

REGISTER_NSP_COMPONENT(ppv) {
    mcc_detail::ss_t nat = par["nat"][0];
    mcc_detail::ss_t pred = par["pred"][0];
    std::cout << mcc_detail::ppv(mcc_detail::get_in(nat, pred)) << std::endl;
}

} // namespace jian
















