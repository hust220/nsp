#ifndef JIAN_CLUSTER_H
#define JIAN_CLUSTER_H

#include "../etl.h"

namespace jian {

class Cluster {
public:
    int _k;
    std::vector<std::vector<int>> _clusters;
    int _max_steps = 200;

    Cluster(int k): _k(k) {}

    void operator ()(const MatrixXf &mat) {
        /// Clear _clusters
        _clusters.clear();
        _clusters.resize(_k);

        /// select k points randomly
        BOOST_ASSERT(mat.rows() == mat.cols());
        rand_select(mat.rows());

    //    for (auto &&cluster: _clusters) {
    //        std::cout << cluster[0] << ' ';
    //    }
    //    std::cout << std::endl;

        int flag = 1;
        int n = _max_steps;
        while (n >= 1 && flag != 0) {
            n--;
            flag = 0;
            /// classify
            for (int i = 0; i < mat.rows(); i++) {
                if (std::any_of(_clusters.begin(), _clusters.end(), [&](const std::vector<int> v){return v[0] == i;}))
                    continue;
                auto min_elem = std::min_element(_clusters.begin(), _clusters.end(), [&](const std::vector<int> &v1, const std::vector<int> &v2){
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
                    return std::accumulate(cluster.begin(), cluster.end(), 0, [&](int sum, int j){
                        return sum + mat(i, j);
                    });
                });
                auto min_elem = std::min_element(cluster.begin(), cluster.end(), [&](int i, int j){
                    return costs[i] < costs[j];
                });
                if (cluster[0] != *min_elem) {
                    flag++;
                    std::swap(cluster[0], *min_elem);
                }
            }
            if (flag == 0 || n == 0) {
    //            for (auto &&cluster: _clusters) {
    //                std::copy(cluster.begin(), cluster.end(), std::ostream_iterator<int>(std::cout, " "));
    //                std::cout << std::endl;
    //            }
            } else {
                for (auto &&cluster: _clusters)
                    cluster.resize(1);
    //            for (auto &&cluster: _clusters) {
    //                std::cout << cluster[0] << ' ';
    //            }
    //            std::cout << std::endl;
            }

        }
    }

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

    void rand_select(int size) {
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

};

} /// namespace jian

#endif




