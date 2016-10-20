#pragma once

#include "Model.hpp"

#include <string>

namespace jian {

class RNA {
public:
    static int rank(const std::string &res_name, const std::string &atom_name);
    static void check(const Model &m);
};

} // namespace jian


