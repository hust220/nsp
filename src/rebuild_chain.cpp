#include <set>
#include <string>
#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom/core.hpp>

namespace jian {

namespace rebuild_chain_detail {

using break_pts_t = std::set<int>;

void set_break_pts(break_pts_t &v, const Chain &c) {
    for (int i = 0; i < c.size() - 1; i++) {
        try {
            if (geom::distance(c[i]["O3*"], c[i+1]["O5*"]) > 6) {
                v.insert(i);
            }    
        } catch(...) {
        }
    }            
}

std::string set_ss(std::string &ss, const break_pts_t &v) {
    std::string str;
    int n = 0;     
    for (auto && c : ss) {
        if (c != '&') {
            str += c;
            if (v.count(n)) {
                str += '&';
            }            
            n++;   
        }            
    }                
    return str;
} // namespace rebuild_chain_detail

}

REGISTER_NSP_COMPONENT(rebuild_chain) {
    rebuild_chain_detail::break_pts_t v;
    rebuild_chain_detail::set_break_pts(v, residues_from_file(par["pdb"][0]));
    std::cout << rebuild_chain_detail::set_ss(par["ss"][0], v);
}

} // namespace jian

