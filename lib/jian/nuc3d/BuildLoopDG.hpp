#pragma once

#include "../utils/Debug.hpp"
#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "../dg.hpp"
#include "C2A.hpp"
#include "HelixPar.hpp"

namespace jian {

class BuildLoopDG {
public:
    Eigen::MatrixXd _dist_bound;
    DihBound _dih_bound;
//    std::deque<std::array<int, 2>> _frags;

    DG dg;

    void init(const std::string &seq, const std::string &ss, const std::string &c = "") {
        int len = seq.size(); _dist_bound.resize(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (i != j) {
                _dist_bound(j, i) = 5; 
                _dist_bound(i, j) = 999;
            } else {
                _dist_bound(i, j) = 0;
            }
        }
        SSTree ss_tree;
        ss_tree.make(seq, ss); 
        LOOP_TRAVERSE(ss_tree.head(), 
            set_bound_loop(_dist_bound, _dih_bound, L); 
            set_bound_helix(_dist_bound, _dih_bound, L->s)
        );
//        set_frags(ss);
    }

//    void set_frags(const std::string &ss) {
//        _frags.clear(); 
//        int beg = 0; 
//        FOR((i, ss.size()-1), IF(ss[i] == '(' && ss[i+1] == ')', _frags.push_back({beg, i}); beg = i+1));
//        _frags.push_back({beg, int(ss.size()-1)});
//    }

    Chain operator ()() {
        Debug::print("## Build Loop By DG\n");
        auto &&c = dg(_dist_bound);
        return c2a(c, 0, c.rows() - 1);
    }

//    Model to_all_atom(const Eigen::MatrixXd &c) {
//        Chain chain;
//        for (auto && i : _frags) {
//            for (auto && j : c2a(c, i[0], i[1])) {
//                chain.push_back(j);
//            }
//        }
//        Model m;
//        m.push_back(chain);
//        return m;
//    }

    void set_bound_loop(Eigen::MatrixXd &b, DihBound &d, loop *l) {
        LOOP_EACH(l, 
            if (RES->next != NULL) {
                if (RES->type == '(' && RES->next->type == ')') {
                    b(RES->num-1, RES->next->num-1) = b(RES->next->num-1, RES->num-1) = HelixPar::dist_bp;
                } else {
                    b(RES->num-1, RES->next->num-1) = b(RES->next->num-1, RES->num-1) = HelixPar::dist_bond;
                }
            }
        );
    }

    void set_bound_helix(Eigen::MatrixXd &b, DihBound &d, const helix &h) {
        int len = 0;
        std::deque<int> s1, s2;
        HELIX_EACH(h, 
            len++;
            s1.push_back(BP->res1.num-1);
            s2.push_back(BP->res2.num-1);
        );
        FOR((i, len), FOR((j, i+1, len), b(s1[i], s1[j]) = b(s1[j], s1[i]) = HelixPar::dist_a(j-i);
                                         b(s2[i], s2[j]) = b(s2[j], s2[i]) = HelixPar::dist_b(j-i);
                                         b(s1[i], s2[j]) = b(s2[j], s1[i]) = HelixPar::dist_c(j-i);
                                         b(s1[j], s2[i]) = b(s2[i], s1[j]) = HelixPar::dist_d(j-i)));
        FOR((i, len), b(s1[i], s2[i]) = b(s2[i], s1[i]) = HelixPar::dist_bp);
    }

};

} // namespace jian

