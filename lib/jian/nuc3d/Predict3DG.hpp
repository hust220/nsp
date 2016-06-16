#pragma once

#include "BasicPredict3D.hpp"
#include "../utils/Par.hpp"

namespace jian {

class Predict3DGImpl;

class Predict3DG : public BasicPredict3D {
public:    
    Predict3DG();
    Predict3DG(const Par &par);
    ~Predict3DG();
    Model predict();
private:
    Predict3DGImpl *_impl;
};

} // namespace jian

