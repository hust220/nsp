#pragma once

#include <memory>

namespace jian {

class Par;
namespace quadruple {class QuaPredImpl;}

class QuaPred {
public:
    QuaPred(const Par &par);
    void predict();

private:
    std::shared_ptr<quadruple::QuaPredImpl> _impl;
};

} // namespace jian


