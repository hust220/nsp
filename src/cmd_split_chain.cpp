#include "nsp.hpp"
#include "pdb.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(split_chain) {
    auto g = par.getv("global");
    auto && m = mol_read_to<Molecule>(g[1]);
    auto nums = par.getv("nums", "n");
    List<Int> v;
    for (auto && n : nums) {
        v.push_back(JN_INT(n));
    }
    for (auto && model : m) {
        Int n = 0;
        Model mod;
        Chain c;
        c.name = "A";
        for (auto && chain : model) {
            for (auto && res : chain) {
                c.push_back(res);
                if (std::find(v.begin(), v.end(), n+1) != v.end()) {
                    if (!c.empty()) {
                        mod.push_back(c);
                        c.clear();
                        c.name = to_str(Char(c.name[0] + 1));
                    }
                }
                n++;
            }
        }
        if (!c.empty()) mod.push_back(c);
        model = mod;
    }
    JN_OUT << m;
}

END_JN

