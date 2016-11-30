#pragma once

#include "Residue.hpp"
#include <iostream>
#include <string>
#include <deque>

namespace jian {

	class Chain : public std::deque<Residue> {
	public:
		std::string name = "A";
		std::string type = "unknown";
		std::string model_name = "unknown";
		std::string m_cg = "aa";

		template<typename _Residues>
		void set_residues(const _Residues &residues) {
			this->clear();
			for (const Residue &res : residues) {
				this->push_back(res);
			}
		}

		JN_DEF_ATOMS;
		JN_DEF_RESIDUES;

	};

#define JN_DEF_CHAINS \
	refs<Chain> chains() { return refs<Chain>().append(*this); }\
	refs<const Chain> chains() const { return refs<const Chain>().append(*this); }


} // namespace jian

