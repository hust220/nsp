#pragma once

#include <string>
#include <deque>
#include "Model.hpp"

namespace jian {

class Molecule : public std::deque<Model> {
public:
    std::string name = "unknown";

    Molecule();
    Molecule(const std::string &file_name);
};

} /// namespace jian

