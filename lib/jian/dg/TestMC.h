#ifndef JIAN_DG_TEST_MC
#define JIAN_DG_TEST_MC

#include "DG.h"

namespace jian {
namespace dg {

class TestMC {
public:
    Log log;

    void operator ()(int n) {
        test_mc(n);
    }

    void test_mc(int n) {
        MatrixXf dist_bound;
        DG::DihBoundType dih_bound;
        std::tie(dist_bound, dih_bound) = bound(n);
        log("test mc:\n");
        DG dg(dist_bound, dih_bound);
        auto c = dg();
        log("coordinate:\n", c, '\n');
    }

    std::pair<DG::DistBoundType, DG::DihBoundType> bound(int n) {
        static auto fn_a = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 2 * 3.14159)) + (2.84 * n) * (2.84 * n));};
        static auto fn_c = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 1.5 * 3.14159)) + (2.84*n-4) * (2.84*n-4));};
        static auto fn_d = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 0.5 * 3.14159)) + (2.84*n+4) * (2.84*n+4));};
        static std::vector<double> chir_data {437, 1498, 2600, 3068, 2442, 696, -1647, -3717, -4580, -3605};

        int len = 2 * (n + 1);
        DG::DistBoundType dist_bound(len, len);
        DG::DihBoundType dih_bound;
//        MatrixXf chir(2 * (n - 2), 5);

//        for (int i = 0; i < len; i++) bound(i, i) = 0;
//
//        for (int i = 0; i < n; i++) {
//            int a = i, b = i + 1, c = 2 * n - i, d = 2 * n + 1 - i;
//            bound(a, d) = bound(d, a) = 15.1;
//            if (i == n - 1) bound(b, c) = bound(c, b) = 15.1;
//        }
//
//        for (int i = 0; i < n - 2; i++) {
//            for (int j = 0; j < 4; j++) {
//                chir(i, j) = i + j;
//                chir(i + n - 2,  j) = i + n + 1 + j;
//            }
//            chir(i, 4) = -11.2116;
//            chir(i + n - 2, 4) = -11.2116;
//        }
//
//        for (int i = 1; i <= n; i++) {
//            for (int j = 0; j <= n - i; j++) {
//                int a = j, b = j + i, c = 2 * n + 1 - j - i, d = 2 * n + 1 - j;
//                bound(a, b) = bound(b, a) = bound(c, d) = bound(d, c) = fn_a(i);
//                bound(a, c) = bound(c, a) = fn_c(i);
//                bound(b, d) = bound(d, b) = fn_d(i);
//            }
//        }

        return {dist_bound, dih_bound};
    }
};



} // namespace dg
} // namespace jian

#endif

