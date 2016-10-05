#pragma once

#include "Dca.hpp"

namespace jian {
namespace dca {

class MfDca : public Dca {
public:
    Matf C, eij;

    void calculate_C();
    virtual void calculate_eij();
    void set_mu(const Matf &m, const Vecf &pi, const Vecf &pj, Vecf &mu1, Vecf &mu2);
    virtual float cal_di(int i, int j);
};

} // namespace dca
} // namespace jian


