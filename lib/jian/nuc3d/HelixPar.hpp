#pragma once

#include <deque>
#include "../utils/traits.hpp"

BEGIN_JN

	struct HelixPar {
		double dih_backbone;
		double dist_bp;
		double dist_bond;
		std::deque<double> dihs;

		HelixPar();
		double dist_a(int n);
		double dist_b(int n);
		double dist_c(int n);
		double dist_d(int n);
		double dih(int n);
	};

END_JN

