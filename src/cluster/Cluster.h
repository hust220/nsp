#ifndef JIAN_CLUSTER_H
#define JIAN_CLUSTER_H

#include <util/util.h>

namespace jian {

class Cluster {
public:
    Cluster(int k): _k(k) {}

    void operator ()(const MatrixXf &mat);
    template<class RandomIt, class Dist>
    void operator ()(RandomIt first, RandomIt last, Dist dist) {
        int size = std::distance(first, last);
        MatrixXf mat(size, size);
        int i, j;
        double distance;
        std::cerr << "Initialization for clustering...       " << std::flush;
        std::cerr << fixed << setprecision(2);
        for (auto it1 = first; it1 != last; it1++) {
            std::cerr << "\b\b\b\b\b\b\b" << setw(6) << std::distance(first, it1) / double(size) * 100.0 << "%" << std::flush;
            for (auto it2 = it1; it2 != last; it2++) {
                i = std::distance(first, it1);
                j = std::distance(first, it2);
                if (i == j) {
                    mat(i, j) = 0;
                } else {
                    distance = dist(*it1, *it2);
                    mat(i, j) = mat(j, i) = distance;
                }
            }
        }
        std::cerr << "\b\b\b\b\b\b\b       \nStart clustering now..." << std::endl;
        (*this)(mat);
    }
    void rand_select(int);

    int _k;
    std::vector<std::vector<int>> _clusters;
    int _max_steps = 200;
};

} /// namespace jian

#endif




