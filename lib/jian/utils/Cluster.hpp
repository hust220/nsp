#pragma once

#include <deque>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>
#include <random>
#include "../matrix.hpp"
#include "log.hpp"
#include "Factory.hpp"
#include "Par.hpp"

#define REG_CLUSTER(name, type) REGISTER_FACTORY(::jian::Cluster::constructor_t, name, type)

BEGIN_JN

	class Cluster {
	public:
		using constructor_t = Cluster *(const Par &);
		using fac_t = Factory<constructor_t>;
		using cluster_t = std::deque<int>;
		using clusters_t = std::deque<cluster_t>;

		clusters_t m_clusters;

		virtual Cluster &operator ()(const Mat & mat) = 0;

		template<typename RandomIt, typename F>
		Cluster &operator ()(RandomIt && first, RandomIt && last, F && dist) {
			Mat *m = to_mat(first, last, dist);
			(*this)(*m);
			delete m;
			return *this;
		}

		template<typename RandomIt, typename F>
		static Mat *to_mat(RandomIt && first, RandomIt && last, F && dist) {
			int i, j, size;
			Num d;
			Mat *mat;

			size = std::distance(first, last);
			mat = new Mat(size, size);
			for (auto it1 = first; it1 != last; it1++) {
				for (auto it2 = it1; it2 != last; it2++) {
					i = std::distance(first, it1);
					j = std::distance(first, it2);
					if (i == j) {
						(*mat)(i, j) = 0;
					}
					else {
						d = dist(*it1, *it2);
						(*mat)(i, j) = (*mat)(j, i) = d;
					}
				}
			}
			return mat;
		}

	};

END_JN

