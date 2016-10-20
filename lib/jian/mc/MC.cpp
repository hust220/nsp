#include "../utils/log.hpp"
#include "MC.hpp"

namespace jian {

MC::MC() {
    Par par(Env::lib() + "/RNA/pars/mc/mc.par");
    par.set(_mc_cycle_steps, "mc_cycle_steps");
    par.set(_mc_write_steps, "mc_write_steps");
    par.set(_mc_heat_steps, "mc_heat_steps");
    par.set(_mc_cool_steps, "mc_cool_steps");
    par.set(_mc_heat_rate, "mc_heat_rate");
    par.set(_mc_init_tempr, "mc_init_tempr");
    par.set(_mc_dec_rate, "mc_dec_rate");
}

bool MC::mc_is_heating() const {
    return _mc_state == 1;
}

bool MC::mc_is_cooling() const {
    return _mc_state == 2;
}

void MC::mc_write() {
    LOG << __FUNCTION__ << "..." << std::endl;
}

void MC::mc_heat() {
    LOG << __FUNCTION__ << "..." << std::endl;
    _mc_state = 1;
    _mc_tempr = _mc_init_tempr;
    mc_base(_mc_heat_steps, [&](double rate){
        if (rate < _mc_heat_rate) {
            _mc_tempr *= 1.1; 
            return true;
        } else if (rate > _mc_heat_rate) {
            _mc_tempr *= 0.9;
            return true;
        } else {
            return true;
        }
    });
}

void MC::mc_cool() {
    LOG << __FUNCTION__ << "..." << std::endl;
    _mc_state = 2;
    int n_rate = 0, n_en = 0;
    double en = 0;
    mc_base(_mc_cool_steps, [&](double rate){
        _mc_tempr *= _mc_dec_rate;
         if (std::fabs(en - _mc_en) <= 0.05) n_en++; else n_en = 0; 
         en = _mc_en;
         if (rate <= 0.1) n_rate++; else n_rate = 0; 
         if (n_rate >= 5 || n_en >= 5) return false; else return true; 
    });
}

void MC::mc() {
    mc_heat();
    mc_cool();
}

double MC::mc_total_energy() {
    return 0;
}

double MC::mc_partial_energy() {
    return rand() * 20;
}

void MC::mc_select() {
    LOG << __FUNCTION__ << "..." << std::endl;
}

void MC::mc_sample() {
    LOG << __FUNCTION__ << "..." << std::endl;
}

void MC::mc_back() {
    LOG << __FUNCTION__ << "..." << std::endl;
}

} // namespace jian

