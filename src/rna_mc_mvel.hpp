#pragma once

namespace jian {

class DHMC;

enum {
    DHMC_MVEL_HELIX,
    DHMC_MVEL_LOOP,
    DHMC_MVEL_FRAG,
    DHMC_MVEL_FRAG3
};

struct dhmc_mvel_t {
    void * p;
    int t;
};

std::ostream &operator <<(std::ostream &stream, const dhmc_mvel_t &el);

void dhmc_mvel_free(dhmc_mvel_t el);

void dhmc_set_mvels(DHMC &dhmc);

void mvels_set_fixed_els(Deque<MvEl *> &mvels, const Deque<Frags> &els);

}

