#include "nsp.hpp"
#include <nsp/nuc2d.hpp>
#include <nsp/dca.hpp>

BEGIN_JN

using namespace dca;

namespace mcc_detail {

// interaction network type
struct in_t {
    double tp, fp, fn;
};

pairs_t get_common_pairs(const pairs_t &pairs1, const pairs_t &pairs2) {
    pairs_t pairs;
    for (auto && pair1 : pairs1) {
        if (std::find_if(pairs2.begin(), pairs2.end(), [&pair1](auto && pair2){
            return pair1[0]==pair2[0]&&pair1[1]==pair2[1];
        }) != pairs2.end()) {
            pairs.push_back(pair1);
        }
    }
    return pairs;
}

/*
pairs_t get_pairs(const ss_t &ss) {
    auto & v = NASS::instance().paired_keys;
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
*/

in_t get_in(const ss_t &nat, const ss_t &pred) {
    in_t in;
    pairs_t pairs_nat = pairs_from_ss(nat);
    pairs_t pairs_pred = pairs_from_ss(pred);
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
    Num n = in.tp + in.fn;
    return n == 0 ? 0 : in.tp / n;
}

double ppv(const in_t &in) {
    Num n = in.tp + in.fp;
    return n == 0 ? 0 : in.tp / n;
}

double mcc(const in_t &in) {
    Num n = sty(in) * ppv(in);
    return n == 0 ? 0 : std::sqrt(n);
}

} // namespace mcc_detail

REGISTER_NSP_COMPONENT(mcc) {
    ss_t nat = par["nat"][0];
    ss_t pred = par["pred"][0];
    std::cout << mcc_detail::mcc(mcc_detail::get_in(nat, pred)) << std::endl;
}

REGISTER_NSP_COMPONENT(sty) {
    ss_t nat = par["nat"][0];
    ss_t pred = par["pred"][0];
    std::cout << mcc_detail::sty(mcc_detail::get_in(nat, pred)) << std::endl;
}

REGISTER_NSP_COMPONENT(ppv) {
    ss_t nat = par["nat"][0];
    ss_t pred = par["pred"][0];
    std::cout << mcc_detail::ppv(mcc_detail::get_in(nat, pred)) << std::endl;
}

REGISTER_NSP_COMPONENT(common_ss) {
    ss_t nat = par["nat"][0];
    ss_t pred = par["pred"][0];
    std::cout << pairs_to_ss(mcc_detail::get_common_pairs(pairs_from_ss(nat), pairs_from_ss(pred)), nat.size()) << std::endl;
}

REGISTER_NSP_COMPONENT(sort_dis) {
    S di_file = par["di"][0];
    int size = std::stoi(par["size"][0]);

    print_tuples(JN_OUT, tuples_from_file(di_file, size));
}

END_JN
















