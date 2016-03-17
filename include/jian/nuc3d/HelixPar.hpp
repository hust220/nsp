#pragma once

#include "../etl.hpp"

namespace jian {

class HelixPar {
public:
    double dih_backbone = 0.264;
    double dist_bp = 15.1;
    double dist_bond = 6.1;
    std::vector<double> _dihs {
        -1.50868, -1.29647, -0.999739, -0.700053, -0.423381, -0.124998, 0.228387, 0.500968,
        0.5183, 0.318411, 0.0464779, -0.18131, -0.314941, -0.345016, -0.285822, -0.163995,
        -0.0100389, 0.138594, 0.234976, 0.243958, 0.168475, 0.0462362, -0.0761236, -0.163417,
        -0.197244, -0.174959, -0.10689, -0.012686, 0.0813646, 0.146262, 0.160388, 0.121046,
        0.0458351, -0.0376064, -0.103726, -0.135776, -0.127817, -0.0843342, -0.0184026, 0.0508864,
        0.102322, 0.119276, 0.0970832, 0.0452772, -0.0174611
    };

    static HelixPar &instance() {
        static HelixPar helix_par;
        return helix_par;
    }

    static double dist_a(int n) {
        return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-2*PI))+(2.84*n)*(2.84*n));
    }

    static double dist_b(int n) {
        return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-2*PI))+(2.84*n)*(2.84*n));
    }

    static double dist_c(int n) {
        return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-1.5*PI))+(2.84*n-4)*(2.84*n-4));
    }

    static double dist_d(int n) {
        return std::sqrt(2*9.7*9.7*(1-std::cos(0.562*n-0.5*PI))+(2.84*n+4)*(2.84*n+4));
    }

    static double dih(int n) {
        return instance()._dihs[n-1];
    }

};

} // namespace jian

