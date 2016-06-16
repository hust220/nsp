#pragma once

#include "BasicPredict3D.hpp"
#include "../utils/Par.hpp"

namespace jian {

class AssembleImpl;

class Assemble : public BasicPredict3D {
public:    
    Assemble();
    Assemble(const Par &par);
    ~Assemble();
    Model predict();
private:
    AssembleImpl *_impl;
};

} // namespace jian

