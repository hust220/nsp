#pragma once

#include <deque>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>
#include <random>

namespace jian {

class Cluster {
public:
    int _k;
    std::deque<std::deque<int>> _clusters;
    int _max_steps = 200;

    Cluster(int k);
    void operator ()(const Eigen::MatrixXd &mat);
    void rand_select(int size);

    template<typename RandomIt, typename F>
    void operator ()(RandomIt first, RandomIt last, F &&dist) {
        int size = std::distance(first, last);
        Eigen::MatrixXd mat(size, size);
        int i, j;
        double distance;
        for (auto it1 = first; it1 != last; it1++) {
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
        (*this)(mat);
    }

};

} // namespace jian

