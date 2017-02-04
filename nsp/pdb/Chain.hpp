#pragma once

#include "Residue.hpp"
#include <iostream>
#include <string>
#include <deque>

BEGIN_JN

	class Chain : public std::deque<Residue> {
	public:
		S name = "A";
		S type = "unknown";
		S model_name = "unknown";
		S m_cg = "aa";

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


END_JN

