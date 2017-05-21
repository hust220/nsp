#pragma once

#include "factory.hpp"
#include "geom.hpp"
#include "cg.hpp"

#define REG_SCORER(name, type) REGISTER_FACTORY(jian::Score::Constructor, name, type)

BEGIN_JN

/**
 * Energy of radius of gyration
 */
template<typename _Chain>
inline Num en_rg_6p(_Chain && c) {
    // Mean position
    Int n = 0;
    Vec v = Vec::Zero(3);
    for (auto && r : c) {
        for (auto && a : r) {
            for (Int i = 0; i < 3; i++) v[i] += a[i];
            n++;
        }
    }
    for (Int i = 0; i < 3; i++) v[i] /= n;

    // sum of distance
    Num e = 0;
    for (auto && r : c) {
        for (auto && a : r) {
            e += geom::dist2(a, v);
        }
    }
    return e / n;
}

/*
template<typename _Chain>
inline Num en_len_6p(const _Chain &c, int beg) {
    Num d, e;

    e = 0;
    if (beg + 1 < c.size()) {
        d = geom::distance(c[beg][0], c[beg][1]);
        e += square(d - m_bond_len_std[1]);
        d = geom::distance(c[beg][1], c[beg + 1][0]);
        e += square(d - m_bond_len_std[1]);
    }
    return e;
}

template<typename _Chain>
inline Num en_ang_6p(const _Chain &c, int beg) {
    Num d, e;

    e = 0;
    if (beg + 2 < c.size()) {
        d = geom::angle(c[beg][0], c[beg + 1][0], c[beg + 2][0]);
        e += square(d - m_bond_angle_std[0]);
    }
    return e;
}

template<typename _Chain>
inline Num en_dih_6p(const _Chain &c, int beg) {
    Num d, e;

    e = 0;
    if (beg + 3 < c.size()) {
        d = geom::dihedral(c[beg][0], c[beg + 1][0], c[beg + 2][0], c[beg + 3][0]);
        e += square(std::sin(0.5*(d - m_bond_dihedral_std[0])));
    }
    if (beg + 1 < c.size()) {
        d = geom::dihedral(c[beg][1], c[beg][0], c[beg + 1][0], c[beg + 1][1]);
        e += square(std::sin(0.5*(d - m_bond_dihedral_std[2])));
    }
    return e;
}

inline Num en_crash_6p(const Residue &r1, const Residue &r2) {
    int i, j;
    Num d, e;
    const Residue *p1, *p2;
    Residue temp1, temp2;

    if (m_cg->is_cg(r1) && m_cg->is_cg(r2)) {
        p1 = &r1;
        p2 = &r2;
    }
    else {
        temp1 = m_cg->to_cg(r1);
        temp2 = m_cg->to_cg(r2);
        p1 = &temp1;
        p2 = &temp2;
    }


    e = 0;
    for (i = 0; i < m_res_size; i++) {
        for (j = 0; j < m_res_size; j++) {
            d = geom::distance(p1->at(i), p2->at(j));
            if (d < 3) {
                e += square(d - 3);
            }
        }
    }

    return e;
}
*/

class Score {
    public:
        using Constructor = Score*(void);
        using fac_t = Factory<Score::Constructor>;

        Num m_score;
        Num m_en_pairing;
        Num m_en_stacking;
        Num m_en_wc;
        Num m_en_nwc;
        Num m_en_vdw;

        CG *m_cg = NULL;

        ~Score() {
            if (m_cg != NULL) delete m_cg;
        }

        virtual void init() = 0;

        virtual void run(const Chain &) = 0;

        virtual void train(const Chain &) = 0;

        virtual void print_counts(std::ostream &) const = 0;

        virtual void print_freqs(std::ostream &) const = 0;

        virtual Num en_stacking(const Residue &r1, const Residue &r2) = 0;

        virtual Num en_pairing(const Residue &r1, const Residue &r2) = 0;

        virtual Score &en_bp(const Residue &r1, const Residue &r2) = 0;

        virtual Num en_len(const Chain &c, int beg) = 0;

        virtual Num en_ang(const Chain &c, int beg) = 0;

        virtual Num en_dih(const Chain &c, int beg) = 0;

        virtual Num en_crash(const Residue &r1, const Residue &r2) = 0;

};

END_JN


