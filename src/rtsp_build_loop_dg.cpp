#include <memory>
#include "cg.hpp"
#include "geom.hpp"
#include "rtsp_build_loop_dg.hpp"

namespace jian {

Chain *build_chain_dg(Str seq, Str ss) {
    static BuildLoopDG builder;
    std::cout << 11 << std::endl;
    builder.init(seq, ss);
    std::cout << 22 << std::endl;
    Chain *chain = new Chain(builder());
    std::cout << 33 << std::endl;
    return chain;
}

void dg_constraints_init(DgConstraints &b, int l) {
    const HelixPar &helix_par = HelixPar::instance();
    for (int i = 0; i < l-1; i++) {
        b.append(DgConstraint{{i, i+1}, {helix_par.dist_bond}});
    }
}

void dg_dist_init(Mat &b, int l) {
    const HelixPar &helix_par = HelixPar::instance();
    b.resize(l, l);
    for (int i = 0; i < l; i++) for (int j = i; j < l; j++) {
        if (j - i == 1) {
            b(i, j) = helix_par.dist_bond;
            b(j, i) = helix_par.dist_bond;
        }
        else if (i != j) {
            b(j, i) = 5;
            b(i, j) = 999;
        }
        else {
            b(i, j) = 0;
        }
    }
}

void dg_constraints_read_loop(DgConstraints &b, const SSE &sse) {
    const HelixPar &helix_par = HelixPar::instance();
    if (sse.has_loop()) {
        auto &loop = sse.loop;
        auto it1 = loop.begin();
        auto it2 = std::next(it1);
        for (; it2 != loop.end(); it1++, it2++) {
            if (it1->type == '(' && it2->type == ')') {
                b.append(DgConstraint{{it1->num - 1, it2->num - 1}, {helix_par.dist_bp}});
            }
            else {
                b.append(DgConstraint{{it1->num - 1, it2->num - 1}, {helix_par.dist_bond}});
            }
        }
    }
}

void dg_constraints_read_helix(DgConstraints &b, const SSE &sse) {
    const HelixPar &helix_par = HelixPar::instance();
    if (sse.has_helix()) {
        auto &helix = sse.helix;
        int len = 0;
        std::deque<int> s1, s2;
        for (auto && bp : helix) {
            len++;
            s1.push_back(bp.res1.num - 1);
            s2.push_back(bp.res2.num - 1);
        }
        for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
            b.append(DgConstraint{{s1[i], s1[j]}, {helix_par.dist_a(j - i)}});
            b.append(DgConstraint{{s2[i], s2[j]}, {helix_par.dist_b(j - i)}});
            b.append(DgConstraint{{s1[i], s2[j]}, {helix_par.dist_c(j - i)}});
            b.append(DgConstraint{{s1[j], s2[i]}, {helix_par.dist_d(j - i)}});
        }
        for (int i = 0; i < len; i++) {
            b.append(DgConstraint{{s1[i], s2[i]}, {helix_par.dist_bp}});
        }
    }
}

void dg_dist_read_loop(Mat &b, const SSE &sse) {
    const HelixPar &helix_par = HelixPar::instance();
    if (sse.has_loop()) {
        auto &loop = sse.loop;
        auto it1 = loop.begin();
        auto it2 = std::next(it1);
        for (; it2 != loop.end(); it1++, it2++) {
            if (it1->type == '(' && it2->type == ')') {
                b(it1->num - 1, it2->num - 1) = b(it2->num - 1, it1->num - 1) = helix_par.dist_bp;
            }
            else {
                b(it1->num - 1, it2->num - 1) = b(it2->num - 1, it1->num - 1) = helix_par.dist_bond;
            }
        }
    }
}

void dg_dist_read_helix(Mat &b, const SSE &sse) {
    const HelixPar &helix_par = HelixPar::instance();
    if (sse.has_helix()) {
        auto &helix = sse.helix;
        int len = 0;
        std::deque<int> s1, s2;
        for (auto && bp : helix) {
            len++;
            s1.push_back(bp.res1.num - 1);
            s2.push_back(bp.res2.num - 1);
        }
        for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
            b(s1[i], s1[j]) = b(s1[j], s1[i]) = helix_par.dist_a(j - i);
            b(s2[i], s2[j]) = b(s2[j], s2[i]) = helix_par.dist_b(j - i);
            b(s1[i], s2[j]) = b(s2[j], s1[i]) = helix_par.dist_c(j - i);
            b(s1[j], s2[i]) = b(s2[i], s1[j]) = helix_par.dist_d(j - i);
        }
        for (int i = 0; i < len; i++) b(s1[i], s2[i]) = b(s2[i], s1[i]) = helix_par.dist_bp;
    }
}

void dg_constraints_read_ss(DgConstraints &constraints, Str seq, Str ss) {
    for (auto &&sse : SSTree(seq, ss, 1)) {
        dg_constraints_read_loop(constraints, sse);
        dg_constraints_read_helix(constraints, sse);
    }
}

void dg_dist_read_ss(Mat &dist, Str seq, Str ss) {
    for (auto &&sse : SSTree(seq, ss, 1)) {
        dg_dist_read_loop(dist, sse);
        dg_dist_read_helix(dist, sse);
    }
}

void dg_constraints_read_chain(DgConstraints &constraints, const Chain &c) {
    JN_TODO;
//    int l = size(c);
//    Vector<Int> v(l);
//    std::iota(v.begin(), v.end(), 0);
//    dg_dist_read_chain(dist, c, v);
}

void dg_dist_read_chain(Mat &dist, const Chain &c) {
    int l = size(c);
    Vector<Int> v(l);
    std::iota(v.begin(), v.end(), 0);
    dg_dist_read_chain(dist, c, v);
}

BuildLoopDG::BuildLoopDG() {
    m_cg.reset(CG::fac_t::create("1p"));
}

BuildLoopDG &BuildLoopDG::init(const Str &seq, const Str &ss) {
    len = size(seq);
    dg_constraints_init(distance_constraints, len);
    dg_constraints_read_ss(distance_constraints, seq, ss);
    return *this;
}

BuildLoopDG &BuildLoopDG::init(const Chain &c, const std::vector<int> &brokens) {
    len = c.size();
    dg_constraints_init(distance_constraints, len);
    dg_constraints_read_brokens(distance_constraints, brokens);
    dg_constraints_read_chain(distance_constraints, c);

    return *this;
}

template<typename T>
static void print_constraints(const T &constraints) {
    for (auto && constraint : constraints) {
        std::cout << join(" ", constraint.key) << ": " << join(" ", constraint.val) << std::endl;
    }
}

Chain BuildLoopDG::operator ()() {
    std::cout << 221 << std::endl;
    print_constraints(distance_constraints);
    auto dg = std::make_unique<Dg>(len, std::array<double, 2>{5, 999}, distance_constraints);
    dg->set_log_file("std.out");
    Mat &&c = dg->sample();
    std::cout << c << std::endl;
    return m_cg->to_aa(c, 0, c.rows() - 1);
}

}

