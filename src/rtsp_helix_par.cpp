#include <vector>
#include <cmath>
#include <string>
#include "rtsp_helix_par.hpp"
#include "math.hpp"
#include "par.hpp"
#include "env.hpp"
#include "jian.hpp"

namespace jian {

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

double HelixPar::dist_a(int n) const {
	return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 2 * PI)) + (2.84*n)*(2.84*n));
}

double HelixPar::dist_b(int n) const {
	return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 2 * PI)) + (2.84*n)*(2.84*n));
}

double HelixPar::dist_c(int n) const {
	return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 1.5*PI)) + (2.84*n - 4)*(2.84*n - 4));
}

double HelixPar::dist_d(int n) const {
	return std::sqrt(2 * 9.7*9.7*(1 - std::cos(0.562*n - 0.5*PI)) + (2.84*n + 4)*(2.84*n + 4));
}

double HelixPar::dih(int n) const {
	return dihs[n - 1];
}

}

