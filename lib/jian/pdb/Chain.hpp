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

void append_chain_to_file(const Chain &chain, const std::string &file_name, int n);
int num_atoms(const Chain &chain);

} // namespace jian

