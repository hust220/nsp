#pragma once

#include <deque>
#include "jian.hpp"

BEGIN_JN

struct HelixPar {
	double dih_backbone;
	double dist_bp;
	double dist_bond;
	std::deque<double> dihs;

	HelixPar();

	static const HelixPar &instance() {
		static HelixPar par;
		return par;
	}
	double dist_a(int n) const;
	double dist_b(int n) const;
	double dist_c(int n) const;
	double dist_d(int n) const;
	double dih(int n) const;
};

END_JN

