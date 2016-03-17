#pragma once

#include "Job.hpp"

namespace jian {

template<typename T>
int triangle_smoothing(T &&b) {
    int f = 0, len = b.rows(); double d;
    FOR((i, len), FOR((j, i+1, len), FOR((k, j+1, len), 
       d = b(i, k) + b(j, k); if (b(i, j) > d) {b(i, j) = d; f++;} if (b(j, i) > d) throw "Wrong bound matrix!";
       d = b(i, j) + b(j, k); if (b(i, k) > d) {b(i, k) = d; f++;} if (b(k, i) > d) throw "Wrong bound matrix!";
       d = b(i, j) + b(i, k); if (b(j, k) > d) {b(j, k) = d; f++;} if (b(k, j) > d) throw "Wrong bound matrix!";
       d = b(k, i) - b(j, k); if (b(j, i) < d) {b(j, i) = d; f++;}
       d = b(k, j) - b(i, k); if (b(j, i) < d) {b(j, i) = d; f++;}
       d = b(j, i) - b(j, k); if (b(k, i) < d) {b(k, i) = d; f++;}
       d = b(k, j) - b(i, j); if (b(k, i) < d) {b(k, i) = d; f++;}
       d = b(j, i) - b(i, k); if (b(k, j) < d) {b(k, j) = d; f++;}
       d = b(k, i) - b(i, j); if (b(k, j) < d) {b(k, j) = d; f++;}
    )));
    return f;
}

template<typename T>
int tetrangle_smoothing(T &&b) {
    int f = 0, len = b.rows(); double d; 
//    FOR((i, len), FOR((j, i+1, len), FOR((k, j+1, len), FOR((l, k+1, len),
//        
//    ))));
    return f;
}

template<typename T>
void smooth(T &&b, double min_dist = 5) {
    // check bound matrix
    if (b.rows() == 0 || b.cols() == 0 || b.rows() != b.cols()) throw "Wrong bound matrix!";

    // check minimum distance
    int len = b.rows();
    FOR((i, b.rows()), FOR((j, b.cols()), IF(i != j && b(i, j) < min_dist, b(i, j) = min_dist)));

    // Main Loop
    int max_step = 100;
    while (triangle_smoothing(b) + tetrangle_smoothing(b) != 0 && max_step-- > 0);
}

} // namespace jian

