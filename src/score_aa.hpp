#pragma once

#include <string>
#include "score_dist.hpp"
#include "score_dih.hpp"
#include "score.hpp"

namespace jian {

class ScoreAa : public Score {
public:
	DistAnal * m_dist_anal = NULL;
	DihAnal  * m_dih_anal = NULL;

	S m_par_dist;
	S m_par_dih;

	double m_bin_dist;
	double m_bin_dih;
	double m_cutoff;

	double m_score_dist;
	double m_score_dih;

	double m_constant;
	double m_weight_dist;
	double m_weight_dih;

	ScoreAa();

	~ScoreAa();

	virtual void init();

	virtual void run(const Chain &);

	virtual void train(const Chain &);

	virtual void print_counts(std::ostream &) const;

	virtual void print_freqs(std::ostream &) const;

	virtual double en_stacking(const Residue &r1, const Residue &r2);

	virtual double en_pairing(const Residue &r1, const Residue &r2);

	virtual Score &en_bp(const Residue &r1, const Residue &r2);

	virtual double en_len(const Chain &c, int beg);

	virtual double en_ang(const Chain &c, int beg);

	virtual double en_dih(const Chain &c, int beg);

	virtual double en_crash(const Residue &r1, const Residue &r2);

};

}

