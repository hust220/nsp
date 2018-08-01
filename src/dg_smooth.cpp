#include "pp.hpp"
#include "matrix.hpp"

namespace jian {

template<typename T>
int dg_triangle_smoothing(T &&b) {
    int f = 0;
    int len = b.rows();
    double d;
    for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) for (int k = j + 1; k < len; k++) {
        d = b(i, k) + b(j, k); if (b(i, j) > d) { b(i, j) = d; f++; } if (b(j, i) > d) throw "Wrong bound matrix!";
        d = b(i, j) + b(j, k); if (b(i, k) > d) { b(i, k) = d; f++; } if (b(k, i) > d) throw "Wrong bound matrix!";
        d = b(i, j) + b(i, k); if (b(j, k) > d) { b(j, k) = d; f++; } if (b(k, j) > d) throw "Wrong bound matrix!";
        d = b(k, i) - b(j, k); if (b(j, i) < d) { b(j, i) = d; f++; }
        d = b(k, j) - b(i, k); if (b(j, i) < d) { b(j, i) = d; f++; }
        d = b(j, i) - b(j, k); if (b(k, i) < d) { b(k, i) = d; f++; }
        d = b(k, j) - b(i, j); if (b(k, i) < d) { b(k, i) = d; f++; }
        d = b(j, i) - b(i, k); if (b(k, j) < d) { b(k, j) = d; f++; }
        d = b(k, i) - b(i, j); if (b(k, j) < d) { b(k, j) = d; f++; }
    }
    return f;
}

template<typename T>
int dg_tetrangle_smoothing(T &&b) {
    // TODO
    return 0;
}

void dg_smooth(Eigen::MatrixXd &b, double min_dist) {
    // check bound matrix
    if (b.rows() == 0 || b.cols() == 0 || b.rows() != b.cols()) throw "Wrong bound matrix!";

    // check minimum distance
    int len = b.rows();
    for (int i = 0; i < b.rows(); i++) for (int j = 0; j < b.cols(); j++) {
        if (i != j && b(i, j) < min_dist) {
            b(i, j) = min_dist;
        }
    }

    // Main Loop
    int max_step = 300;
    while (dg_triangle_smoothing(b) + dg_tetrangle_smoothing(b) != 0 && max_step-- > 0);
}

}

