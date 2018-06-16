#pragma once

#include <string>
#include <deque>
#include "pdb_model.hpp"

namespace jian {

class Molecule : public std::deque<Model> {
public:
    S name = "unknown";
	S m_cg = "aa";

	JN_DEF_ATOMS;
	JN_DEF_RESIDUES;
	JN_DEF_CHAINS;
	JN_DEF_MODELS;

};


} /// namespace jian

