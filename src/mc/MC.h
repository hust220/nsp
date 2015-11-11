#ifndef MC_H
#define MC_H

#include <util/util.h>
#include "System.h"

namespace jian {

namespace mc {
    
class MC {
public:
    MC(): unif_distr(0, 1) {}
    void operator ()(System &);

    std::random_device rd;
    std::uniform_real_distribution<> unif_distr;

};

} /// namespace mc

} /// namespace jian


#endif




