#pragma once

#include <string>
#include <deque>
#include "Model.hpp"

namespace jian {

class Molecule : public std::deque<Model> {
public:
    std::string name = "unknown";
	std::string m_cg = "aa";

	JN_DEF_ATOMS;
	JN_DEF_RESIDUES;
	JN_DEF_CHAINS;
	JN_DEF_MODELS;

};


} /// namespace jian

