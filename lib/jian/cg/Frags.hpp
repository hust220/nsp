#pragma once

#include <deque>
#include <string>
#include <map>

#include "../pdb.hpp"

namespace jian {

class CG;

	class Frags {
	public:
		CG *m_cg;
		std::vector<Mat *> m_mats;
		std::vector<Chain *> m_chains_aa;
		std::vector<Chain *> m_chains_cg;
		std::map<std::string, std::vector<int>> m_ids;
		std::string m_path;
		int m_frag_size;
		int m_res_size;
		static std::map<std::string, std::map<int, Frags>> m_instances;

		static Frags &instance(std::string cg, int frag_size);

		~Frags();

		void init(std::string cg, int frag_size);

		Chain get_chain(int i, const Mat &c);

		//template<typename Coord, typename Frag>
		//auto run(Coord &&coord, Frag &&frag) {
		//	int i, j, k, num_atoms, len, index;
		//	Chain residues;

		//	num_atoms = m_frag_size * m_res_size;
		//	len = frag[1] - frag[0] + 1;
		//	Chain chain;
		//	for (i = 0; i < len - num_atoms + 1; i += m_res_size) {
		//		Mat c(num_atoms, 3);
		//		for (j = 0; j < num_atoms; j++) {
		//			for (k = 0; k < 3; k++) {
		//				c(j, k) = coord(frag[0] + i + j, k);
		//			}
		//		}
		//		std::deque<double> scores;
		//		for (auto && frag : m_frags) {
		//			auto r = geom::suppos(*frag, c);
		//			scores.push_back(r.rmsd);
		//		}
		//		auto min = std::min_element(scores.begin(), scores.end());
		//		index = std::distance(scores.begin(), min);
		//		residues = get_chain(index, c);
		//		if (i == 0) {
		//			for (j = 0; j < m_frag_size - 1; j++) {
		//				chain.push_back(residues[j]);
		//			}
		//		}
		//		chain.push_back(residues[m_frag_size - 1]);
		//	}
		//	return chain;
		//}

		void extract_frags(const std::string &pdb);
	};
}