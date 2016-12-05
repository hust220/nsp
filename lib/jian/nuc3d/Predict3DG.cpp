#include "Predict3DG.hpp"
#include "Predict3DGImpl.hpp"

BEGIN_JN

Predict3DG::Predict3DG() : _impl(NULL) {}

Predict3DG::Predict3DG(const Par &par) : _impl(new Predict3DGImpl(par)) {}

Predict3DG::~Predict3DG() {
    if (_impl != NULL) delete _impl;
}

Model Predict3DG::predict() {
    return _impl->predict();
}

}

