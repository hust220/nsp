#include "DBSCAN.hpp"

BEGIN_JN
	REG_CLUSTER("dbscan", DBSCAN);

	DBSCAN::DBSCAN(const Par &par) {
		m_eps = 1.0;
		m_min_pts = 10;

		par.set(m_eps, "eps");
		par.set(m_min_pts, "min_pts");
	}

	void DBSCAN::init(int l, const Mat &mat) {
		m_is_visited.resize(l);
		for (auto && b : m_is_visited) b = false;
		m_neighbors.resize(l);
		for (int i = 0; i < l; i++) {
			for (int j = i + 1; j < l; j++) {
				if (mat(i, j) < m_eps) {
					m_neighbors[i].push_back(j);
					m_neighbors[j].push_back(i);
				}
			}
		}
		//LOG << mat << std::endl;
		//LOG << m_eps << ' ' << m_min_pts << std::endl;
		//for (auto && neighbors : m_neighbors) {
		//	for (auto && i : neighbors) LOG << i << ' ';
		//	LOG << std::endl;
		//}
		m_is_noise.resize(l);
		for (auto && b : m_is_noise) b = false;
	}

	void DBSCAN::add_to_cluster(int i, cluster_t &cluster) {
		cluster.push_back(i);
		for (auto && j : m_neighbors[i]) {
			if (!m_is_visited[j]) {
				m_is_visited[j] = true;
				if (m_neighbors[j].size() >= m_min_pts) {
					add_to_cluster(j, cluster);
				}
				else {
					m_is_noise[j] = true;
				}
			}
		}
	}

	Cluster &DBSCAN::operator ()(const Mat & mat) {
		int i, l;

		if (mat.rows() != mat.cols()) {
			throw "jian::DBSCAN::operator() error! The rows and the cols should be equal!";
		}

		l = mat.rows();
		init(l, mat);
		for (i = 0; i < l; i++) {
			if (m_is_visited[i]) {
				continue;
			}
			else {
				m_is_visited[i] = true;
				if (m_neighbors[i].size() < m_min_pts) {
					m_is_noise[i] = true;
				}
				else {
					cluster_t cluster;
					add_to_cluster(i, cluster);
					m_clusters.push_back(std::move(cluster));
				}
			}
		}
		for (i = 0; i < l; i++) {
			if (m_is_noise[i]) {
				//cluster_t cluster;
				//cluster.push_back(i);
				m_clusters.push_back({ i });
			}
		}

		std::sort(m_clusters.begin(), m_clusters.end(), [](auto &&c1, auto &&c2) {return c1.size() > c2.size(); });

		return *this;
	}


}
