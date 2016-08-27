#include "BuildLoopDG.hpp"

namespace jian {

Chain *build_chain_dg(std::string seq, std::string ss) {
    static BuildLoopDG builder;
    builder.init(seq, ss);
    Chain *chain = new Chain(std::move(builder()));
    return chain;
}

void BuildLoopDG::init(const std::string &seq, const std::string &ss, const std::string &c) {
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
}

Chain BuildLoopDG::operator ()() {
    Debug::print("## Build Loop By DG\n");
    auto &&c = dg(_dist_bound);
    return c2a(c, 0, c.rows() - 1);
}

void BuildLoopDG::set_bound_loop(Eigen::MatrixXd &b, DihBound &d, loop *l) {
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

void BuildLoopDG::set_bound_helix(Eigen::MatrixXd &b, DihBound &d, const helix &h) {
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

} // namespace jian

