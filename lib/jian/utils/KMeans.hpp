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

	class KMeans : public Cluster {
	public:
		int m_k;
		int m_max_steps;

		KMeans(const Par &par);

		void rand_select(int size);

		virtual Cluster &operator ()(const Mat &mat);

	};

END_JN

