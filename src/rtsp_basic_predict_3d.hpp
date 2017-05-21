#pragma once

#include "pdb_model.hpp"

BEGIN_JN

class BasicPredict3D {
public:
    virtual Model predict() = 0;
};

END_JN

