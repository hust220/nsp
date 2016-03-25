#include "nsp.hpp"
#include <jian/nuc3d/Assemble.hpp>
#include <jian/nuc3d/Predict3FA.hpp>
#include <jian/nuc3d/Predict3DG.hpp>
#include <jian/nuc3d/Predict3MC.hpp>
#include <jian/nuc3d/triass/TriAss.hpp>

namespace jian {

template<typename T, typename U>
double residue_distance(T &&r1, U &&r2) {
    return geom::distance(r1["C4*"], r2["C4*"]);
}

template<typename T, typename U>
bool satisfy(T &&model, U &&c) {
    double cutoff = 20;
    EACH_RES1(model, EACH_RES2(model, 
        if (has_contact(c, N_RES1, N_RES2) && residue_distance(RES1, RES2) > cutoff) return false));
    return true;
}

REGISTER_NSP_COMPONENT(3drna) {
    static std::map<std::string, std::function<BasicPredict3D*(const Par &)>> _3drna_methods = {
        {"FA", [](const Par &par){return new Assemble(par);}},
        {"FADG", [](const Par &par){return new Assemble(par);}},
        {"DG", [](const Par &par){return new Predict3DG(par);}},
        {"tri-assemble", [](const Par &par){return new TriAss(par);}}
    };

    auto job = JobPredict3D(par);
    auto c = read_constraints(job._file_constraints);
    display_start_information_job(job);
    auto method = _3drna_methods[job._method](par);
    int max_it_num = 100; int step = 1; int num = 0;
    while (true) {
        auto &&model = method->predict();
        Trace::log("Model ", step, " : ");
        if (satisfy(model, c)) {
            write_pdb(model, job._name + "-" + JSTR(num+1) + ".pdb");
            Trace::log("Satisfy\n");
            num++;
        } else {
            Trace::log("Unsatisfy\n");
        }
        if (num >= job._num || step >= max_it_num) {
            break;
        }
        step++;
    }
    delete method;
    display_end_information_job(job);
}

REGISTER_NSP_COMPONENT(r3d) {
    Predict3MC r3d(par);
    r3d.predict();
}

REGISTER_NSP_COMPONENT(mc) {
    MC mc;
    mc.mc_heat();
    mc.mc_cool();
}

} // namespace jian

