#pragma once

#include "NucSS.hpp"
#include "../pdb/IFModel.hpp"

namespace jian {
namespace nuc2d {

class GetSS {
public:
    typedef std::pair<int, int> BP;
    typedef std::list<BP> BPs;

    template<typename ModelType> 
    std::string operator ()(ModelType &&model) {
        return get_ss(model);
    }

    template<typename ModelType> 
    std::string get_ss(ModelType &&model) {
        auto bps = get_bps(model);
//        print_bps(bps);
        int len = pdb::num_residues(model);
        return bps_to_ss(bps, len);
    }

    void print_bps(const BPs &bps) const {
        for (auto && bp : bps) {
            std::cout << bp.first << '-' << bp.second << std::endl;
        }
    }

    template<typename MolType> 
    BPs get_bps(MolType &&model) {
        BPs bps;
        int num_res1 = 0;
        for (auto &&chain : model) for (auto &&residue : chain) {
            int num_res2 = 0;
            for (auto &&chain2 : model) for (auto &&residue2 : chain2) {
                if (num_res2 > num_res1) {
                    if (pdb::is_residue_paired(residue, residue2)) {
                        append(bps, std::make_pair(num_res1, num_res2));
                    }
                }
                num_res2++;
            }
            num_res1++;
        }
        return bps;
    }

    std::string bps_to_ss(const BPs &bps, int len) {
        std::string ss(len, '.');
        for (auto && bp : bps) {
            auto level = bp_level(ss, bp);
            if (level >= 0) {
                ss[bp.first] = NucSS::instance().paired_keys[level].first;
                ss[bp.second] = NucSS::instance().paired_keys[level].second;
            }
        }
        return ss;
    }

    int bp_level(const std::string &ss, const BP &bp) {
        int left = bp.first, right = bp.second;
        if (ss[left] != '.' or ss[right] != '.') return -1;
        for (int i = 0; i < NucSS::instance().paired_keys.size(); i++) {
            int score = 0, flag = 1;
            for (int j = left + 1; j < right; j++) {
                if (ss[j] == NucSS::instance().paired_keys[i].first) score++;
                else if (ss[j] == NucSS::instance().paired_keys[i].second) score--;
                if (score < 0) break;
            }
            if (score != 0) flag = 0;
            if (flag == 1) return i;
        }
        return -1;
    }

};

} // namespace nuc2d
} // namespace jian

