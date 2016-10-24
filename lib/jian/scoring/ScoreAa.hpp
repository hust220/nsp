#pragma once

#include "Score.hpp"
#include "DistAnal.hpp"
#include "DihAnal.hpp"
#include "../cg.hpp"

namespace jian {

	struct AA {};

	template<>
	class Score<AA> : public ScoreBase {
	public:
		DistAnal * m_dist_anal = NULL;
		DihAnal  * m_dih_anal = NULL;

		std::string m_par_dist;
		std::string m_par_dih;

		double m_bin_dist = 0.3;
		double m_bin_dih = 4.5;
		double m_cutoff = 20;

		double m_score_dist = 0;
		double m_score_dih = 0;

		double m_constant = 27.1118;
		double m_weight_dist = 0.433513;
		double m_weight_dih = 1.59348;


		~Score();
		virtual void init();
		virtual void run(const Chain &);
		virtual void train(const Chain &);
		virtual void print_counts(std::ostream &) const;
	};

	using ScoreAa = Score<AA>;

} // namespace jian

