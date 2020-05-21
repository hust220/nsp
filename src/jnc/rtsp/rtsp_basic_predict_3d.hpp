#pragma once

#include "pdb_model.hpp"

namespace jian {

class BasicPredict3D {
public:
    virtual Model predict() = 0;
};

}

