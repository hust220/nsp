#include <cassert>
#include <algorithm>
#include <numeric>
#include "log.hpp"
#include "serial.hpp"
#include "mc.hpp"

namespace jian {

    MC::MC() {
        Par par(Env::lib() + "/RNA/pars/mc/mc.par");
        par.set(_mc_cycle_steps, "mc_cycle_steps");
        par.set(_mc_write_steps, "mc_write_steps");
        par.set(_mc_heat_steps, "mc_heat_steps");
        par.set(_mc_cool_steps, "mc_cool_steps");
        par.set(_mc_heat_rate, "mc_heat_rate");
        _mc_init_tempr = 20;
        par.set(_mc_init_tempr, "mc_init_tempr");
        _mc_inc_rate = 1.05;
        par.set(_mc_inc_rate, "mc_inc_rate");
        _mc_dec_rate = 0.9995;
        par.set(_mc_dec_rate, "mc_dec_rate");
        _mc_lowest_en = 0.5;
        par.set(_mc_lowest_en, "mc_lowest_en");
        _mc_lowest_rate = 0.01;
        par.set(_mc_lowest_rate, "mc_lowest_rate");
        _mc_lowest_tempr = 20;
        par.set(_mc_lowest_tempr, "mc_lowest_tempr");
        _mc_highest_tempr = 200;
        par.set(_mc_highest_tempr, "mc_highest_tempr");
    }

    void MC::mc_next_step() {
        _mc_step++;
    }

    bool MC::mc_is_heating() const {
        return _mc_state == MC_HEATING;
    }

    bool MC::mc_is_cooling() const {
        return _mc_state == MC_COOLING;
    }

    void MC::mc_write() {
        LOG << __FUNCTION__ << "..." << std::endl;
    }

    void MC::mc_run() {
        tokenize_v v, w, u;
        int steps;

        //LOG << __FUNCTION__ << "..." << std::endl;
        tokenize(_mc_queue, v, "+");
        for (auto && s : v) {
            tokenize(s, w, ":");
            if (w.size() >= 2) {
                if (w[0] == "heat") {
                    steps = JN_INT(w[1]);
                    if (w.size() == 3) _mc_init_tempr = JN_DBL(w[2]);
                    mc_heat(steps);
                }
                else if (w[0] == "cool") {
                    steps = JN_INT(w[1]);
                    if (w.size() == 3) _mc_init_tempr = JN_DBL(w[2]);
                    mc_cool(steps);
                }
                else if (w[0] == "warm") {
                    steps = JN_INT(w[1]);
                    if (w.size() == 3) _mc_init_tempr = JN_DBL(w[2]);
                    mc_warm(steps);
                }
                else if (w[0] == "samc") {
                    steps = JN_INT(w[1]);
                    if (w.size() >= 3) {
                        tokenize(w[2], u, "-");
                        if (size(u) >= 2) {
                            _mc_highest_tempr = JN_DBL(u[0]);
                            _mc_lowest_tempr = JN_DBL(u[1]);
                        }
                        else if (size(u) == 1) {
                            _mc_highest_tempr = JN_DBL(u[0]);
                        }
                        else {
                            throw "Illegal settings of SAMC temperature!";
                        }
                    }
                    if (w.size() >= 4) {
                        _mc_dec_rate = JN_NUM(w[3]);
                    }
                    _mc_init_tempr = _mc_highest_tempr;
                    _mc_num_samc = 1;
                    if (size(w) >= 4) _mc_num_samc = JN_INT(w[3]);
                    for (int i = 0; i < _mc_num_samc; i++) {
                        mc_samc(steps);
                    }
                }
                else if (w[0] == "remc") {
                    steps = JN_INT(w[1]);
                    tokenize(w[2], u, "-");
                    if (u.size() != 2) {
                        throw "Illegal settings of REMC temperature!";
                    }
                    _mc_lowest_tempr = JN_DBL(u[0]);
                    _mc_highest_tempr = JN_DBL(u[1]);
                    mc_remc(steps);
                }
                else {
                    throw w[0] + ": illegal mc mode!";
                }
            }
            else {
                throw "jian::MC::mc_run error!";
            }
        }
    }

    void MC::mc_heat(int steps) {
        //LOG << __FUNCTION__ << "..." << std::endl;
        _mc_tempr = _mc_init_tempr;
        _mc_state = MC_HEATING;
        mc_base(steps, [this]() {
            if (_mc_local_succ_rate < _mc_heat_rate) {
                _mc_tempr *= 1.1;
                return true;
            }
            else if (_mc_local_succ_rate > _mc_heat_rate) {
                _mc_tempr *= 0.9;
                return true;
            }
            else {
                return true;
            }
        });
    }

    void MC::mc_cool(int steps) {
        //LOG << __FUNCTION__ << "..." << std::endl;
        _mc_tempr = _mc_init_tempr;
        _mc_state = MC_COOLING;
        int n_rate = 0, n_en = 0;
        double en = 0;
        mc_base(steps, [this, &en, &n_en, &n_rate]() {
            _mc_tempr *= _mc_dec_rate;
            if (std::fabs(en - _mc_en) <= 0.05) n_en++; else n_en = 0;
            en = _mc_en;
            if (_mc_local_succ_rate <= 0.05) n_rate++; else n_rate = 0;
            return n_rate < 5 && n_en < 5;
        });
    }

    void MC::mc_warm(int steps) {
        //LOG << __FUNCTION__ << "..." << std::endl;
        _mc_tempr = _mc_init_tempr;
        _mc_state = MC_WARMING;
        mc_base(steps, [this]() {
            return true;
        });
    }

    void MC::mc_samc(int steps) {
        //LOG << __FUNCTION__ << "..." << std::endl;
        _mc_tempr = _mc_init_tempr;
        _mc_state = MC_SAMC;
        int n_rate = 0, n_en = 0;
        double en = 0;
        mc_base(steps, [this, &en, &n_en, &n_rate]() {
            _mc_tempr *= _mc_dec_rate;
            n_en += (std::fabs(en - _mc_en) <= _mc_lowest_en ? 1 : 0);
            en = _mc_en;
            n_rate += (_mc_local_succ_rate <= _mc_lowest_rate ? 1 : 0);
            return n_rate < 10 && n_en < 10 && _mc_tempr >= _mc_lowest_tempr;
        });
    }

    void MC::mc_remc(int steps) {
        //LOG << __FUNCTION__ << "..." << std::endl;
#ifndef JN_PARA
        throw "REMC is only supported in the parallel version! "
            "You may need to recompile nsp if you really want to use REMC!";
#else
        _mc_state = MC_REMC;
        if (mpi_size() == 1) {
            _mc_tempr = _mc_lowest_tempr;
        }
        else {
            _mc_tempr = _mc_lowest_tempr + mpi_rank() * (_mc_highest_tempr - _mc_lowest_tempr) / (mpi_size() - 1);
        }
        Serial serial;

        auto exchange_tempr = [](auto &&tempr, auto &&en) {
            int l = tempr.size();
            double d;

            for (int i = 0; i < l - 1; i++) {
                d = en[i] - en[i + 1];
                d *= 1.0 / tempr[i] - 1.0 / tempr[i + 1];
                d = std::exp(d);
                if (rand() < std::min(1.0, d)) {
                    std::swap(tempr[i], tempr[i+1]);
                }
            }
        };

        mc_base(steps, [this, &serial, &exchange_tempr]() {
            if (mpi_rank() == 0) {
                std::vector<double> tempr(mpi_size());
                std::vector<double> en(mpi_size());
                tempr[0] = _mc_tempr;
                en[0] = _mc_en;
                for (int i = 1; i < mpi_size(); i++) {
                    serial.parse(mpi_recv(i), tempr[i], en[i]);
                }
                exchange_tempr(tempr, en);
                _mc_tempr = tempr[0];
                for (int i = 1; i < mpi_size(); i++) {
                    mpi_send(serial.stringify(tempr[i]), i);
                }
            }
            else {
                mpi_send(serial.stringify(_mc_tempr, _mc_en), 0);
                serial.parse(mpi_recv(0), _mc_tempr);
            }
            return true;
        });
#endif
    }

//    double MC::mc_total_energy() {
//        return 0;
//    }
//
//    double MC::mc_partial_energy() {
//        return rand() * 20;
//    }
//
//    void MC::mc_select() {
//        //LOG << __FUNCTION__ << "..." << std::endl;
//    }
//
//    void MC::mc_sample() {
//        //LOG << __FUNCTION__ << "..." << std::endl;
//    }
//
//    void MC::mc_backup() {
//        //LOG << __FUNCTION__ << "..." << std::endl;
//    }
//
//    void MC::mc_rollback() {
//        //LOG << __FUNCTION__ << "..." << std::endl;
//    }

} // namespace jian

