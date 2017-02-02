#pragma once

#include "Model.hpp"

#include <string>

BEGIN_JN

class RNA {
public:
    static int rank(const S &res_name, const S &atom_name);
    static void check(const Model &m);
};

END_JN


