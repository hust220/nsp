#pragma once

#include "Score.hpp"

BEGIN_JN
	class Score6p : public Score {
	public:
		Score6p();

		virtual double en_len(const Chain &c, int beg);

		virtual double en_ang(const Chain &c, int beg);

		virtual double en_dih(const Chain &c, int beg);

		virtual double en_crash(const Residue &r1, const Residue &r2);

		virtual bool in_base(int);
	};
}