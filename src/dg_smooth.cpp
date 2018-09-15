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

bool dg_smooth(Eigen::MatrixXd &b, std::array<int, 2> &pair) {
    // check bound matrix
    if (b.rows() == 0 || b.cols() == 0 || b.rows() != b.cols()) throw "Wrong bound matrix!";

    // Main Loop
    int max_step = 300;
    while (dg_triangle_smoothing(b) + dg_tetrangle_smoothing(b) != 0 && max_step-- > 0);

    // Find out the pair that has the maximum distance
    double max = -1;
    double d;
    int len = b.rows();
    for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
        d = b(i, j) - b(j, i);
        if (d > max) {
            max = d;
            pair[0] = i;
            pair[1] = j;
        }
    }

    return max == 0;
}

}

