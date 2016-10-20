#pragma once

#include "Residue.hpp"
#include <iostream>
#include <string>
#include <deque>

namespace jian {

class Chain : public std::deque<Residue> {
public:
    std::string name = "X";
    std::string type = "unknown";
    std::string model_name = "unknown";
};

} // namespace jian

