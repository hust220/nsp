#pragma once

namespace jian {

struct HelixPar {
    static double dih_backbone;
    static double dist_bp;
    static double dist_bond;
    static double dist_a(int n);
    static double dist_b(int n);
    static double dist_c(int n);
    static double dist_d(int n);
    static double dih(int n);
};

} // namespace jian

