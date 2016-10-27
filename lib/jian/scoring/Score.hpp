#pragma once

#include "../utils/Factory.hpp"

#define REG_SCORER(name, type) REGISTER_FACTORY(jian::ScoreBase::Constructor, name, type)

namespace jian {

	class ScoreBase {
	public:
		using Constructor = ScoreBase *(void);
		//typedef ScoreBase *(*Constructor)();

		double m_score;

		virtual void init() = 0;

		virtual void run(const Chain &) = 0;

		virtual void train(const Chain &) = 0;

		virtual void print_counts(std::ostream &) const = 0;

		virtual void print_freqs(std::ostream &) const = 0;

		virtual double en_stacking(const Residue &r1, const Residue &r2) = 0;

		virtual double en_pairing(const Residue &r1, const Residue &r2) = 0;

	};

	using FacScorer = Factory<ScoreBase::Constructor>;

	template<typename CG_T>
	class Score;

}

