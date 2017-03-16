#pragma once

#include <jian/utils/Factory.hpp>
#include <jian/geom.hpp>
#include "../cg.hpp"

#define REG_SCORER(name, type) REGISTER_FACTORY(jian::Score::Constructor, name, type)

BEGIN_JN

/**
 * Energy of radius of gyration
 */
template<typename _Chain>
Num en_rg(_Chain && c) {
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

class Score {
public:
	using Constructor = Score*(void);
	using fac_t = Factory<Score::Constructor>;

	double m_score;
	double m_en_pairing;
	double m_en_stacking;
	double m_en_wc;
	double m_en_nwc;
	double m_en_vdw;

	CG *m_cg = NULL;

	~Score() {
		if (m_cg != NULL) delete m_cg;
	}

	virtual void init() = 0;

	virtual void run(const Chain &) = 0;

	virtual void train(const Chain &) = 0;

	virtual void print_counts(std::ostream &) const = 0;

	virtual void print_freqs(std::ostream &) const = 0;

	virtual double en_stacking(const Residue &r1, const Residue &r2) = 0;

	virtual double en_pairing(const Residue &r1, const Residue &r2) = 0;

	virtual Score &en_bp(const Residue &r1, const Residue &r2) = 0;

	virtual double en_len(const Chain &c, int beg) = 0;

	virtual double en_ang(const Chain &c, int beg) = 0;

	virtual double en_dih(const Chain &c, int beg) = 0;

	virtual double en_crash(const Residue &r1, const Residue &r2) = 0;

};

END_JN


