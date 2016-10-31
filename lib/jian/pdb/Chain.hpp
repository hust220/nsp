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
};

} // namespace jian

