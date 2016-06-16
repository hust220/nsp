#include "MC.hpp"

namespace jian {

MC::MC() {
    Par par(Env::lib() + "/RNA/pars/mc/mc.par");
    par.set(_mc_cycle_steps, "mc_cycle_steps");
    par.set(_mc_write_steps, "mc_write_steps");
    par.set(_mc_heat_steps, "mc_heat_steps");
    par.set(_mc_cool_steps, "mc_cool_steps");
    par.set(_mc_heat_rate, "mc_heat_rate");
    par.set(_mc_tempr, "mc_init_tempr");
}

bool MC::mc_is_heating() const {
    return _mc_state == 1;
}

bool MC::mc_is_cooling() const {
    return _mc_state == 2;
}

void MC::mc_write() {
    std::cout << __FUNCTION__ << "..." << std::endl;
}

void MC::mc_heat() {
    std::cout << __FUNCTION__ << "..." << std::endl;
    _mc_state = 1;
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
    std::cout << __FUNCTION__ << "..." << std::endl;
    _mc_state = 2;
    int temp = 0, flag = 0;
    double en = 0;
    mc_base(_mc_cool_steps, [&](double rate){
        _mc_tempr *= 0.999;
        if (std::fabs(en - _mc_en) <= 0.05) {
            flag++;
        } else {
            flag = 0;
        }
        en = _mc_en;
        if (rate <= 0.01) {
            temp++;
        } else {
            temp = 0;
        }
        if (temp >= 5 || flag >= 5) {
            return false;
        } else {
            return true;
        }
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
    std::cout << __FUNCTION__ << "..." << std::endl;
}

void MC::mc_sample() {
    std::cout << __FUNCTION__ << "..." << std::endl;
}

void MC::mc_back() {
    std::cout << __FUNCTION__ << "..." << std::endl;
}

} // namespace jian

