#pragma once

#include "../utils/Factory.hpp"
#include "../cg.hpp"

#define REG_SCORER(name, type) REGISTER_FACTORY(jian::ScoreBase::Constructor, name, type)

namespace jian {

	class ScoreBase {
	public:
		using Constructor = ScoreBase*(void);
		using fac_t = Factory<ScoreBase::Constructor>;

		double m_score;
		double m_en_pairing;
		double m_en_stacking;
		double m_en_wc;
		double m_en_nwc;
		double m_en_vdw;

		CG *m_cg = NULL;

		~ScoreBase() {
			if (m_cg != NULL) delete m_cg;
		}

		virtual void init() = 0;

		virtual void run(const Chain &) = 0;

		virtual void train(const Chain &) = 0;

		virtual void print_counts(std::ostream &) const = 0;

		virtual void print_freqs(std::ostream &) const = 0;

		virtual double en_stacking(const Residue &r1, const Residue &r2) = 0;

		virtual double en_pairing(const Residue &r1, const Residue &r2) = 0;

		virtual ScoreBase &en_bp(const Residue &r1, const Residue &r2) = 0;

		virtual double en_len(const Chain &c, int beg) = 0;

		virtual double en_ang(const Chain &c, int beg) = 0;

		virtual double en_dih(const Chain &c, int beg) = 0;

		virtual double en_crash(const Residue &r1, const Residue &r2) = 0;

	};

}

