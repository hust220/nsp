#pragma once

#include <deque>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>
#include <random>
#include "../utils/log.hpp"

namespace jian {

class Cluster {
public:
    using result_t = std::deque<std::deque<int>>;

    int _k;
    result_t _clusters;
    int _max_steps = 200;

    Cluster(int k);
    result_t operator ()(const Eigen::MatrixXd &mat);
    void rand_select(int size);

    template<typename RandomIt, typename F>
    result_t operator ()(RandomIt first, RandomIt last, F &&dist) {
        int size = std::distance(first, last);
		LOG << "size: " << size << std::endl;
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
        return (*this)(mat);
    }

};

} // namespace jian

