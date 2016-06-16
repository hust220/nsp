#include <set>
#include <string>
#include <iostream>
#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/geom/core.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(rebuild_chain) {
    auto rna = RNA(par["pdb"][0]);
    std::set<int> break_points;
    Residue old_res; 
    int res_num = 0; 
    for (auto &&chain: rna) {
        for (auto &&residue: chain) {
            if (! old_res.empty()) {
                try {
                    if (geom::distance(residue["O5*"], old_res["O3*"]) > 3.5) {
                        break_points.insert(res_num - 1);
                    }    
                } catch(...) {
                }
            }        
            old_res = residue;
            res_num++;
        }            
    }                
    std::string ss;  
    res_num = 0;     
    for (int i = 0; i < par["ss"][0].size(); i++) {
        if (par["ss"][0][i] == '&') {
            continue;
        }            
        ss += par["ss"][0][i];
        if (break_points.count(res_num)) {
            ss += '&';
        }            
        res_num++;   
    }                
    std::cout << ss;
}

} // namespace jian

