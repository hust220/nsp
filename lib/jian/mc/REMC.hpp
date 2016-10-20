#pragma once

#include <cmath>
#include <iostream>
#include "../utils/rand.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"

namespace jian {

void remc_init(double min, double max, int n);

template<typename T>
void remc_run(T &&sys) {
    sys.run();
}

class REMC {
public:
    double _mc_tempr = 20;
    int _mc_cycle_steps;
    int _mc_write_steps;
    int _mc_step = 0;
    double _mc_local_succ_rate;
    double _mc_en;
    double _mc_heat_rate;
    int _mc_state = 0; // 0: ready, 1: heating, 2: cooling, 3: done

    REMC();

    void run();
    virtual void mc_write();
    virtual double mc_total_energy();
    virtual double mc_partial_energy();
    virtual void mc_select();
    virtual void mc_sample();
    virtual void mc_back();
    bool mc_is_heating() const;
    bool mc_is_cooling() const;

};

} // namespace jian

