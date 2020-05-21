#pragma once

#include "rss_sst.hpp"

namespace jian {

namespace nuc2d {

class MergeRings {
public:
    void operator ()(SST &sst) {
        merge_rings(sst);
    }

    void merge_rings(SST &sst) {
        S seq, ss;
        std::tie(seq, ss) = get_seq_ss(sst);
        std::cout << seq << "\n" << ss << std::endl;
        auto ssts = get_constraint_trees(sst, seq, ss);
//        for (auto &&temp_sst: ssts) {
//            std::cout << "SST:\n";
//            temp_sst.apply([](const SST &temp_sst){std::cout << temp_sst._data << std::endl;});    
//            std::cout << std::endl;
//        }
        for (auto &&temp_sst: ssts) temp_sst.apply([&](const SST &temp_sst){merge_ring(sst, temp_sst._data);});
    }

    std::pair<std::string, std::string> get_seq_ss(const SST &sst) {
        S seq, ss;
        int size = sst.fold([](const SSE &hairpin){return hairpin.size();}, std::plus<int>());
        (seq.reserve(size), ss.reserve(size));
        seq += sst._data._helix[0].seq() + sst._data._loop.front().seq();
        ss += sst._data._helix[0].ss() + sst._data._loop.front().ss();
        for (int i = 0; i < sst._sons.size(); i++) {
            S temp_seq, temp_ss;
            std::tie(temp_seq, temp_ss) = get_seq_ss(*std::next(sst._sons.begin(), i));
            seq += temp_seq + std::next(sst._data._loop.begin(), i + 1)->seq();
            ss += temp_ss + std::next(sst._data._loop.begin(), i + 1)->ss();
        }
        seq += sst._data._helix[1].seq();
        ss += sst._data._helix[1].ss();
        return std::make_pair(seq, ss);
    }

    std::list<SST> get_constraint_trees(const SST &sst, const S seq, const S ss) {
        std::list<SST> list;
        BuildSST build_sst;
        for (int i = 1; i < NASS::instance().paired_keys.size(); i++) {
            char left = NASS::instance().paired_keys[i].first, right = NASS::instance().paired_keys[i].second;
            if (!std::count_if(ss.begin(), ss.end(), [&](const char &c){return c == left || c == right;})) break;
            build_sst.set_primary_keys(left, right);
            list.push_back(build_sst(seq, ss));
        }
        return list;
    }

    void merge_ring(SST &sst, const SSE &hairpin) {
        if (hairpin._helix[0].empty()) return;
        SST *left, *right;
        std::tie(left, right) = find_left_right(sst, hairpin);
        std::cout << left->_data << right->_data << std::endl;
        auto ring = find_ring(sst, hairpin);
        auto pseudo_knot = make_pseudo_knot(ring);
    }

    std::pair<SST *, SST *> find_left_right(SST &sst, const SSE &hairpin) {
        std::pair<SST *, SST *> pair;
        sst.apply([&](SST &temp_sst){
            for (auto &&frag: temp_sst._data._loop) {
                if (has_same_residues(frag, hairpin._helix[0])) pair.first = temp_sst.self();
                if (has_same_residues(frag, hairpin._helix[1])) pair.second = temp_sst.self();
            }
        });
        return pair;
    }

    bool has_same_residues(const Frag &frag1, const Frag &frag2) {
        for (auto &&res1: frag1) for (auto &&res2: frag2) if (res1 == res2) return true;
        return false;
    }

    Tree<SST *> find_ring(const SST &sst, const SSE &hairpin) {}

    SST make_pseudo_knot(const Tree<SST *> &rint) {
        return SST();
    }
};

} // namespace nuc2d
} // jian

