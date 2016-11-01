#pragma once

#include <string>
#include "DistAnal.hpp"
#include "DihAnal.hpp"
#include "ScoreBase.hpp"

namespace jian {

	class ScoreAa : public ScoreBase {
	public:
		DistAnal * m_dist_anal = NULL;
		DihAnal  * m_dih_anal = NULL;

		std::string m_par_dist;
		std::string m_par_dih;

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

		virtual ScoreBase &en_bp(const Residue &r1, const Residue &r2);

		virtual double en_len(const Residue &r1, const Residue &r2);

		virtual double en_ang(const Residue &r1, const Residue &r2, const Residue &r3);

		virtual double en_dih(const Residue &r1, const Residue &r2, const Residue &r3, const Residue &r4);

		virtual double en_crash(const Residue &r1, const Residue &r2);

	};

} // namespace jian

