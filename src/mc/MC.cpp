#include "MC.h"

namespace jian {

namespace mc {

void MC::operator ()(System &sys) {
    double E = sys.energy();
    int succ_num = 0;
    double succ_rate = 1;
    int total_steps = 10000;
    int echo_steps = 1000;
    int cycle_steps = 50;
    double tempr = 20;
    int local_succ_num = 0;
    double en_atom_old = 0;
    double en_atom_new = 0;

    /* print initial condition */
    std::cerr << "=============== Simulated Annealing =================" << std::endl;
    std::cerr << "Initial condition." << std::endl;
    std::cerr << "Energy: " << sys.energy() << std::endl;
    std::cerr << "Temprature: " << tempr << std::endl;

    /* warming */
    std::cerr << "Step 1: warming... " << std::endl;
    std::cerr << setw(8) << "step" << setw(15) << "energy" << setw(15) << "temprature" << setw(15) << "success rate" << std::endl;
    int warming_steps = 1000;
    for (int i = 0;; i++) {
        en_atom_old = sys.energy();
        sys.move();
        en_atom_new = sys.energy();
        if (en_atom_new > en_atom_old) { // randomly reject or accept
            if (unif_distr(rd) > exp(-(en_atom_new - en_atom_old) / tempr)) { // reject
                sys.rollback();
            } else { // accept
                E = E - en_atom_old + en_atom_new;
                local_succ_num++;
            }
        } else { // accept
            E = E - en_atom_old + en_atom_new;
            local_succ_num++;
        }

        if (i % cycle_steps == cycle_steps - 1) {
            double local_succ_rate = double(local_succ_num) / cycle_steps;
            if (local_succ_rate < 0.8) {
                tempr *= 2;
            } else {
                std::cerr << setw(8) << i + 1 << setw(15) << E << setw(15) << tempr << setw(15) << local_succ_rate << std::endl;
                break;
            }
            local_succ_num = 0;
            //succ_rate = succ_num / (i + 1.0);
            std::cerr << setw(8) << i + 1 << setw(15) << E << setw(15) << tempr << setw(15) << local_succ_rate << std::endl;
        }
    }

    /* cooling */
    std::cerr << "Step 2: cooling... " << std::endl;
    std::cerr << setw(8) << "step" << setw(15) << "energy" << setw(15) << "temprature" << setw(15) << "success rate" << std::endl;
    for (int i = 0; i < total_steps; i++) {
        en_atom_old = sys.energy();
        sys.move();
        en_atom_new = sys.energy();
//        int m = int(unif_distr(rd) * len);
//        int n = int(unif_distr(rd) * 3);
//        double x = (unif_distr(rd) - 0.5) * 2;
//        en_atom_old = atom_energy(c, m);
//        c(m, n) += x;
//        en_atom_new = atom_energy(c, m);
        if (en_atom_new > en_atom_old) { // randomly reject or accept
            if (unif_distr(rd) > exp(-(en_atom_new - en_atom_old) / tempr)) { // reject
//                c(m, n) -= x;
                sys.rollback();
            } else { // accept
                E = E - en_atom_old + en_atom_new;
                local_succ_num++;
            }
        } else { // accept
            E = E - en_atom_old + en_atom_new;
            local_succ_num++;
            
            /// store minimal energy and correspondingly coordinates
            if (E < sys.min_energy()) {
                sys.set_min_state();
            }
        }

        if (i % cycle_steps == cycle_steps - 1) {
            double local_succ_rate = double(local_succ_num) / cycle_steps;
            //if (local_succ_rate > 0.5) {
            //    tempr *= 0.5;
            //} else {
                tempr *= 0.9;
            //}
            local_succ_num = 0;
            //succ_rate = succ_num / (i + 1.0);
            std::cerr << setw(8) << i + 1 << setw(15) << E << setw(15) << tempr << setw(15) << local_succ_rate << std::endl;
        }
    }

    std::cerr << "minimal E: " << sys.min_energy() << ' ' << E << std::endl;
}

} /// namespace mc

} /// namespace jian


