#pragma once

#include "JobPredict3D.hpp"
#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "../dg.hpp"
#include "C2A.hpp"
#include "HelixPar.hpp"
#include "../utils/Debug.hpp"
#include "../utils/rand.hpp"
#include "../utils/file.hpp"

namespace jian {

class Predict3DGImpl : public BasicPredict3D, public JobPredict3D {
public:
    Eigen::MatrixXd _dist_bound;
    DihBound _dih_bound;
    std::deque<std::array<int, 2>> _frags;
    double _min_dist = 5;
    double _max_dist = 999;

    DG dg;

    Predict3DGImpl(const Par &par) : JobPredict3D(par) {
//        par.set(dg._num_mc_steps, "n_mc");
        Debug::print("# Set Bound Eigen::Matrix\n");
        set_bound();
        set_frags();
    }

    void set_bound() {
        int len = _seq.size();
        init_bound(len);
        set_bound_constraints();
        for (auto && pair : NucSS::instance().paired_keys) {
            auto ss = _ss;
            for (auto && c : ss) {
                if (c == pair.first) {
                    c = '(';
                } else if (c == pair.second) {
                    c = ')';
                } else {
                    c = '.';
                }
            }
            Debug::print(ss, '\n');
            if (std::any_of(ss.begin(), ss.end(), [](auto &&c){return c != '.';})) {
                SSTree ss_tree; 
                ss_tree.make(_seq, ss);
                LOOP_TRAVERSE(ss_tree.head(), 
                    set_bound_loop(_dist_bound, _dih_bound, L); 
                    set_bound_helix(_dist_bound, _dih_bound, L->s)
                );
            } else {
                break;
            }
        }
    }

    void init_bound(int len) {
        _dist_bound.resize(len, len); 
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (i != j) {
                _dist_bound(j, i) = _min_dist; 
                _dist_bound(i, j) = _max_dist;
            } else {
                _dist_bound(i, j) = 0;
            }
        }
    }

    void set_frags() {
        int beg = 0; 
        for (int i = 0; i < _ss.size() - 1; i++) {
            if (_ss[i] == '(' && _ss[i+1] == ')') {
                _frags.push_back({beg, i}); 
                beg = i + 1;
            }
        }
        _frags.push_back({beg, int(_ss.size()-1)});
    }

    void set_bound_constraints() {
        EACH_SPLIT_LINE(_file_constraints.c_str(), " ",
            if (F.size() == 3) {
                int i = JN_INT(F[0]) - 1;
                int j = JN_INT(F[1]) - 1;
                double d = JN_DBL(F[2]);
                _dist_bound(i, j) = std::min(_max_dist, d + 8);
                _dist_bound(j, i) = std::max(_min_dist, d - 8);
            }
        );
    }

    Model predict() {
        Debug::print("# Predict 3D Structure by DG\n");
        Debug::print(_dist_bound, '\n');
        auto &&c = dg(_dist_bound, _dih_bound);
        Debug::print(c, '\n');
        return to_all_atom(c);
    }

    Model to_all_atom(const Eigen::MatrixXd &c) {
        Chain chain; EACH(i, _frags, EACH(j, c2a(c, i[0], i[1]), chain.push_back(j)));
        Model m; m.push_back(chain); 
        return m;
    }

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
            s2.push_back(BP->res2.num-1)
        );
        for (int i = 0; i < len; i++) for (int j = 0; j < i + 1; j++) {
            b(s1[i], s1[j]) = b(s1[j], s1[i]) = HelixPar::dist_a(j-i);
            b(s2[i], s2[j]) = b(s2[j], s2[i]) = HelixPar::dist_b(j-i);
            b(s1[i], s2[j]) = b(s2[j], s1[i]) = HelixPar::dist_c(j-i);
            b(s1[j], s2[i]) = b(s2[i], s1[j]) = HelixPar::dist_d(j-i);
            d[{s1[i], s1[j], s2[j], s2[i]}] = HelixPar::dih(j-i);
        }
        FOR((i, len), b(s1[i], s2[i]) = b(s2[i], s1[i]) = HelixPar::dist_bp);
    }

};

} // namespace jian
