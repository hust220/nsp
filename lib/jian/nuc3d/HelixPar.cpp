#include <vector>
#include <cmath>
#include <string>
#include "HelixPar.hpp"
#include "../utils/math.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"

namespace jian {

class HelixParImpl {
public:
    double dih_backbone;
    double dist_bp;
    double dist_bond;
    std::deque<double> dihs;

    HelixParImpl() {
        std::string path = Env::lib() + "/RNA/pars/nuc3d/HelixPar/helix.pars";
        Par par(path);
        dih_backbone = std::stod(par["dih_backbone"][0]);
        dist_bp = std::stod(par["dist_bp"][0]);
        dist_bond = std::stod(par["dist_bond"][0]);
        for (auto && s : par["dihs"]) {
            dihs.push_back(std::stod(s));
        }
    }
};

HelixParImpl l_helix_par_impl;

double HelixPar::dih_backbone = l_helix_par_impl.dih_backbone;

double HelixPar::dist_bp = l_helix_par_impl.dist_bp;

double HelixPar::dist_bond = l_helix_par_impl.dist_bond;

double HelixPar::dist_a(int n) {
    return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-2*PI))+(2.84*n)*(2.84*n));
}

double HelixPar::dist_b(int n) {
    return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-2*PI))+(2.84*n)*(2.84*n));
}

double HelixPar::dist_c(int n) {
    return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-1.5*PI))+(2.84*n-4)*(2.84*n-4));
}

double HelixPar::dist_d(int n) {
    return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-0.5*PI))+(2.84*n+4)*(2.84*n+4));
}

double HelixPar::dih(int n) {
    return l_helix_par_impl.dihs[n-1];
}

} // namespace jian

