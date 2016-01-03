#ifndef JIAN_MC_SYSTEM_H
#define JIAN_MC_SYSTEM_H

namespace jian {

namespace mc {

class System {
public:    
    virtual void move() = 0;
    virtual void rollback() = 0;
    virtual double energy() = 0;
    virtual void set_min_state() = 0;
    virtual double min_energy() = 0;
};

} /// namespace mc

} /// namespace jian

#endif

