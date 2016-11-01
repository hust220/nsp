#pragma once

#include "Score.hpp"

namespace jian {
	class Score6p : public Score {
	public:
		Score6p();

		virtual double en_len(const Residue &r1, const Residue &r2);

		virtual double en_ang(const Residue &r1, const Residue &r2, const Residue &r3);

		virtual double en_dih(const Residue &r1, const Residue &r2, const Residue &r3, const Residue &r4);

		virtual double en_crash(const Residue &r1, const Residue &r2);

		virtual bool in_base(int);
	};
}