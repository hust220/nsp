#include "KMeans.hpp"

BEGIN_JN
	REG_CLUSTER("kmeans", KMeans);

	KMeans::KMeans(const Par &par) {
		m_k = 1;
		m_max_steps = 200;

		par.set(m_k, "k");
		par.set(m_max_steps, "max_steps");
	}

	void KMeans::rand_select(int size) {
		std::vector<int> a(size);
		std::iota(a.begin(), a.end(), 0);

		std::random_device rd;
		std::mt19937 gen(rd());

		for (int i = 0; i < m_k; i++) {
			std::uniform_int_distribution<> dis(0, a.size() - 1);
			int n = dis(gen);
			m_clusters[i].push_back(a[n]);
			a[n] = a.back();
			a.pop_back();
		}
	}

	Cluster &KMeans::operator ()(const Mat &mat) {
		// Clear m_clusters
		m_clusters.clear();
		m_clusters.resize(m_k);

		// select k points randomly
		rand_select(mat.rows());

		int flag = 1;
		int n = m_max_steps;
		while (n >= 1 && flag != 0) {
			n--;
			flag = 0;
			/// classify
			for (int i = 0; i < mat.rows(); i++) {
				if (std::any_of(m_clusters.begin(), m_clusters.end(), [&](auto &&v) {return v[0] == i; }))
					continue;
				auto min_elem = std::min_element(m_clusters.begin(), m_clusters.end(), [&](auto &&v1, auto &&v2) {
					Num d1 = mat(i, v1[0]);
					Num d2 = mat(i, v2[0]);
					return d1 < d2;
				});
				min_elem->push_back(i);
			}

			/// select new center
			for (auto &&cluster : m_clusters) {
				std::vector<Num> costs(cluster.size());
				std::transform(cluster.begin(), cluster.end(), costs.begin(), [&](int i) {
					return std::accumulate(cluster.begin(), cluster.end(), 0.0, [&](Num sum, int j)->Num {
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
			}
			else {
				for (auto &&cluster : m_clusters)
					cluster.resize(1);
			}
		}

		std::sort(m_clusters.begin(), m_clusters.end(), [](auto &&c1, auto &&c2) {return c1.size() > c2.size(); });

		return *this;
	}



}
