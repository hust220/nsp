#pragma once

#include <memory>

namespace jian {

class Par;
namespace triple {class TriPredImpl;}

class TriPred {
public:
    TriPred(const Par &par);
    void predict();

private:
    std::shared_ptr<triple::TriPredImpl> _impl;
};

} // namespace jian


