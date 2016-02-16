#ifndef JIAN_NUC3D_HELIX_PAR
#define JIAN_NUC3D_HELIX_PAR

#include "../etl.h"

namespace jian {
namespace nuc3d {

class HelixPar {
public:
    static constexpr double dih_backbone = 0.264;

    static double dist_a(int n) {
        return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 2 * 3.14159)) + (2.84 * n) * (2.84 * n));
    }

    static double dist_c(int n) {
        return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 1.5 * 3.14159)) + (2.84*n-4) * (2.84*n-4));
    }

    static double dist_d(int n) {
        return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 0.5 * 3.14159)) + (2.84*n+4) * (2.84*n+4));
    }
};

} // namespace nuc3d
} // namespace jian

#endif

