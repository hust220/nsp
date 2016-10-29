#pragma once

#include <string>
#include <deque>
#include "Model.hpp"

namespace jian {

class Molecule : public std::deque<Model> {
public:
    std::string name = "unknown";

	template<typename CG_T>
	bool is_cg() const {
		return m_cg == CG_T::m_cg;
	}

	template<typename CG_T>
	Molecule & cg() {
		for (auto && model : *this) {
			model.cg<CG_T>();
		}
		m_cg = CG_T::m_cg;
		return *this;
	}

private:
	Residue::cg_code m_cg = Residue::CG_AA;

};


} /// namespace jian

