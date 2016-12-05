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
#include "Cluster.hpp"

BEGIN_JN

	class DBSCAN : public Cluster {
	public:
		Num m_eps;
		int m_min_pts;
		std::vector<int> m_is_visited;
		std::vector<std::list<int>> m_neighbors;
		std::vector<bool> m_is_noise;

		DBSCAN(const Par &par);

		void init(int l, const Mat &mat);

		void add_to_cluster(int i, cluster_t &cluster);

		virtual Cluster &operator ()(const Mat & mat);

	};

END_JN

