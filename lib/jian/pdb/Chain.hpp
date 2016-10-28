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

	template<typename CG_T>
	bool is_cg() const {
		return m_cg == CG_T::m_cg;
	}

	template<typename CG_T>
	Chain & cg() {
		*this = CG_T::chain(*this);
		m_cg = CG_T::m_cg;
		return *this;
	}

private:
	Residue::cg_code m_cg = Residue::CG_AA;
};

} // namespace jian

