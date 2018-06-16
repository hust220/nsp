#pragma once

#include <deque>
#include <string>
#include <map>

#include "pdb.hpp"

namespace jian {

class CG;

class ResFrags {
public:
	CG *m_cg;
	std::vector<Mat *> m_mats;
	std::vector<Chain *> m_chains_aa;
	std::vector<Chain *> m_chains_cg;
	std::map<std::string, std::vector<int>> m_ids;
	S m_path;
	int m_frag_size;
	int m_res_size;
	static std::map<std::string, std::map<int, ResFrags>> m_instances;

	static ResFrags &instance(S cg, int frag_size);

	~ResFrags();

	void init(S cg, int frag_size);

	Chain get_chain(int i, const Mat &c);

	void extract_frags(const S &pdb);
};
}
