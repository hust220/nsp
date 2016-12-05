#include <vector>
#include <cmath>
#include <string>
#include "HelixPar.hpp"
#include "../utils/math.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"
#include "../utils/traits.hpp"

BEGIN_JN


	HelixPar::HelixPar() {
		S path = Env::lib() + "/RNA/pars/nuc3d/HelixPar/helix.pars";
		Par par(path);
		dih_backbone = std::stod(par["dih_backbone"][0]);
		dist_bp = std::stod(par["dist_bp"][0]);
		dist_bond = std::stod(par["dist_bond"][0]);
		for (auto && s : par["dihs"]) {
			dihs.push_back(std::stod(s));
		}
	}

	double HelixPar::dist_a(int n) {
		return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 2 * PI)) + (2.84*n)*(2.84*n));
	}

	double HelixPar::dist_b(int n) {
		return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 2 * PI)) + (2.84*n)*(2.84*n));
	}

	double HelixPar::dist_c(int n) {
		return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 1.5*PI)) + (2.84*n - 4)*(2.84*n - 4));
	}

	double HelixPar::dist_d(int n) {
		return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 0.5*PI)) + (2.84*n + 4)*(2.84*n + 4));
	}

	double HelixPar::dih(int n) {
		return dihs[n - 1];
	}

END_JN

