#ifndef JIAN_DG_MC_H
#define JIAN_DG_MC_H

#include "Gradient.hpp"
#include "Move.hpp"

namespace jian {
namespace dg {

// simulated annealing monte carlo method

class MC : public virtual Gradient {
public:
    double _mc_tempr = 20;

    MC() = default;
    MC(const MC &) = default;
    MC &operator =(const MC &) = default;

    template<typename Fn> 
    void base_mc(int steps, Fn &&ctrl_tempr) {
        Move mv;
        int cycle_steps = 10, len = c.rows(), local_succ_num = 0;
        double en = total_energy(c), min_en = en;
        auto min_coord = c;
        for (int i = 0; i < steps; i++) {
            auto index = mv.pick(len);
            auto en_atom_old = atom_energy(c, index);
            mv.move(c);
            auto en_atom_new = atom_energy(c, index);
            if (en_atom_new > en_atom_old and rand() > exp(-(en_atom_new - en_atom_old) / _mc_tempr)) {
                mv.back(c);
            } else {
                en = en - en_atom_old + en_atom_new;
                local_succ_num++;
                if (en < min_en) (min_en = en, min_coord = c);
            }

            if (i % cycle_steps == cycle_steps - 1) {
                double local_succ_rate = double(local_succ_num) / cycle_steps;
                local_succ_num = 0;
                log(i + 1, ": ", total_dist_energy(c), "(dist energy) ", total_dih_energy(), "(dih energy) ", 
                    _mc_tempr, "(temprature) ", local_succ_rate, "(success rate)\n");
                if (not ctrl_tempr(local_succ_rate)) break;
            }
        }
        c = min_coord;
        gradient();
    }

    void heat() {
        log("Step 1: heating...\n");
        base_mc(50, [&](double rate){if (rate < 0.8) {_mc_tempr *= 2; return true;} else return false;});
    }

    void cool() {
        log("Step 2: cooling...\n");
        base_mc(_num_mc_steps, [&](double rate){
            if (rate > 0.5) {
                _mc_tempr *= 0.9; 
            } else if (rate > 0.2) {
                _mc_tempr *= 0.999;
            } else if (rate > 0.005) {
                _mc_tempr *= 0.99999;
            }
            return true;
        });
    }

    void mc() {
        init_dihs();
        log("Start MC...\n",
            "Energy: ", total_dist_energy(c), "(dist) ", total_dih_energy(), "(dih)\n",
            "Temprature: ", _mc_tempr, '\n');

        heat();
        cool();

        log("Finish MC.\nDistance energy: ", total_dist_energy(c), " Chirality energy: ", total_dih_energy(), '\n');
    }

};

} // namespace dg
} // namespace jian

#endif


