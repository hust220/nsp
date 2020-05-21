#pragma once

#include "pdb_model.hpp"

#include <string>

namespace jian {

class RNA {
public:
    static int rank(const S &res_name, const S &atom_name);
    static void check(const Model &m);
};

}


