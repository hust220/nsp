#include <deque>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>
#include <random>
#include "Cluster.hpp"

namespace jian {

Cluster::Cluster(int k): _k(k) {}

Cluster::result_t Cluster::operator ()(const Eigen::MatrixXd &mat) {
    // Clear _clusters
    _clusters.clear();
    _clusters.resize(_k);

    // select k points randomly
    rand_select(mat.rows());

    int flag = 1;
    int n = _max_steps;
    while (n >= 1 && flag != 0) {
        n--;
        flag = 0;
        /// classify
        for (int i = 0; i < mat.rows(); i++) {
            if (std::any_of(_clusters.begin(), _clusters.end(), [&](auto &&v){return v[0] == i;}))
                continue;
            auto min_elem = std::min_element(_clusters.begin(), _clusters.end(), [&](auto &&v1, auto &&v2){
                double d1 = mat(i, v1[0]);
                double d2 = mat(i, v2[0]);
                return d1 < d2;
            });
            min_elem->push_back(i);
        }

        /// select new center
        for (auto &&cluster: _clusters) {
            std::vector<double> costs(cluster.size());
            std::transform(cluster.begin(), cluster.end(), costs.begin(), [&](int i){
                return std::accumulate(cluster.begin(), cluster.end(), 0.0, [&](double sum, int j)->double{
                    return sum + mat(i, j);
                });
            });
            auto min_elem = std::min_element(costs.begin(), costs.end());
            if (min_elem != costs.begin()) {
                flag++;
                std::swap(cluster[0], cluster[std::distance(costs.begin(), min_elem)]);
            }
        }
        if (flag == 0 || n == 0) {
        } else {
            for (auto &&cluster: _clusters)
                cluster.resize(1);
        }

    }

    std::sort(_clusters.begin(), _clusters.end(), [](auto &&c1, auto &&c2){return c1.size() > c2.size();});

    return _clusters;
}

void Cluster::rand_select(int size) {
    std::vector<int> a(size);
    std::iota(a.begin(), a.end(), 0);

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < _k; i++) {
        std::uniform_int_distribution<> dis(0, a.size() - 1);
        int n = dis(gen);
        _clusters[i].push_back(a[n]);
        a[n] = a.back();
        a.pop_back();
    }
}

} // namespace jian

