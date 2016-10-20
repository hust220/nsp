#pragma once

#include "../pdb/Model.hpp"

namespace jian {

class BasicPredict3D {
public:
    virtual Model predict() = 0;
};

} // namespace jian

