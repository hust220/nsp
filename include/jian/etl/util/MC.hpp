#pragma once

#include "../std.hpp"
#include "Rand.hpp"

namespace jian {

class MC : public virtual Rand {
public:
    double _mc_tempr = 20;
    int _mc_cycle_steps = 20;
    int _mc_write_steps = 100;
    int _mc_step = 0;
    double _mc_local_succ_rate;
    double _mc_en;
    int _mc_heat_steps = 100;
    int _mc_cool_steps = 10000;

    template<typename Fn> 
    void base_mc(int steps, Fn &&ctrl_tempr) {
        int local_succ_num = 0;
        _mc_en = mc_energy();
        mc_write();
        for (_mc_step = 0; _mc_step < steps; _mc_step++) {
            mc_select();
            auto &&en_old = mc_partial_energy();
            mc_sample();
            auto &&en_new = mc_partial_energy();
            auto &&en_diff = en_new - en_old;
            if (en_new > en_old && rand() > std::exp(-en_diff / _mc_tempr)) {
                mc_back();
            } else {
                _mc_en += en_diff;
                local_succ_num++;
            }

            if (_mc_step % _mc_write_steps == _mc_write_steps - 1) {
                mc_write();
            }
            if (_mc_step % _mc_cycle_steps == _mc_cycle_steps - 1) {
                _mc_local_succ_rate = double(local_succ_num) / _mc_cycle_steps;
                local_succ_num = 0;
                if (! ctrl_tempr(_mc_local_succ_rate)) break;
            }
        }
    }

    virtual void mc_write() {
        std::cout << __FUNCTION__ << "..." << std::endl;
    }

    virtual void mc_heat() {
        std::cout << __FUNCTION__ << "..." << std::endl;
        base_mc(_mc_heat_steps, [&](double rate){
            if (rate < 0.5) {
                _mc_tempr *= 2; 
            }
            return true;
        });
    }

    virtual void mc_cool() {
        std::cout << __FUNCTION__ << "..." << std::endl;
        base_mc(_mc_cool_steps, [&](double rate){
            if (rate > 0.5) {
                _mc_tempr *= 0.99; 
            } else if (rate > 0.2) {
                _mc_tempr *= 0.999;
            } else {
                _mc_tempr *= 0.9999;
            }
            return true;
        });
    }

    virtual double mc_energy() {
        return rand() * 20;
    }

    virtual double mc_partial_energy() {
        return rand() * 20;
    }

    virtual void mc_select() {
        std::cout << __FUNCTION__ << "..." << std::endl;
    }

    virtual void mc_sample() {
        std::cout << __FUNCTION__ << "..." << std::endl;
    }

    virtual void mc_back() {
        std::cout << __FUNCTION__ << "..." << std::endl;
    }

};

} // namespace jian

