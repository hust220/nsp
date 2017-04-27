#include "bps.hpp"
#include "bires.hpp"
#include "bires_models.hpp"
#include "../scoring/ParBp.hpp"
#include "../nuc2d/get_ss.hpp"
#include <jian/utils/Env.hpp>
#include <jian/geom.hpp>

BEGIN_JN

static Bool is_close(const Residue &r1, const Residue &r2) {
    return geom::distance(r1["O3*"], r2["C5*"]) < 4;
}

/*
static Bps model_bps(const Model &m) {
    Bps bps;
    auto && rs = m.residues();
    Int l = size(rs);
    for (Int i = 0; i < l; i++) {
        for (Int j = i + 1; j < l; j++) {
            if (is_bp(rs[i], rs[j])) {
                bps.push_back({i, j});
            }
        }
    }
    return bps;
}

static Str model_ss(const Model &m) {
    auto && bps = model_bps(m);
    Int l = num_residues(m);
    return bps_to_ss(bps, l);
}
*/

struct BiResModels {

    Array<Model, 16*4> models;

    BiResModels() {
        for_each_model(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"), [this](const Model &m, int n) {
                extract(m);
                });
        Int i = 0;
        /*
        for (auto && c1 : {'A', 'U', 'G', 'C'}) {
            Int j = 0;
            for (auto && c2 : {'A', 'U', 'G', 'C'}) {
                std::cout << c1 << c2 << ' ';
                for (Int k = 0; k < 4; k++) {
                    std::cout << size(models[(i*4+j)*4+k]) << ' ';
                }
                std::cout << std::endl;
                j++;
            }
            i++;
        }
        */
    }

    void extract(const Model &m) {
        Chain dq;
        Deque<Char> dq_ss;

        Str seq = JN_ seq(m);
        //std::cout << seq << std::endl;
        Str ss = get_ss(m.residues());
        //std::cout << ss << std::endl;
        auto bps = ss_to_bps(ss);
        auto rs = m.residues();
        for (auto && bp : bps) {
            dq.push_back(rs[bp[0]]);
            dq.push_back(rs[bp[1]]);
            //std::cout << bp[0] << ' ' << bp[1] << ' ' << geom::distance(rs[bp[0]]["C1*"], rs[bp[1]]["C1*"]) << std::endl;
            models[bires_code(seq[bp[0]], seq[bp[1]], ss[bp[0]], ss[bp[1]])].push_back(dq);
            dq.clear();
        }
        //std::cout << ss << std::endl;
        Int l = size(ss);
        for (Int i = 0; i < l-1; i++) {
            if (is_close(rs[i], rs[i+1])) {
                dq.push_back(rs[i]);
                dq.push_back(rs[i+1]);
                models[bires_code(seq[i], seq[i+1], ss[i], ss[i+1])].push_back(dq);
                dq.clear();
            }
        }
        /*
           Int n = 0;
           for (const Residue & r : m.residues()) {
           if (!dq.empty() && !is_close(dq.back(), r)) {
           dq.clear();
           }
           dq.push_back(r);
           if (size(dq) == 2) {
           std::cout << n-1 << ' ' << n << ' ' << geom::distance(rs[n-1]["C1*"], rs[n]["C1*"]) << std::endl;
           models[bires_code(seq[n-1], seq[n], ss[n-1], ss[n])].push_back(dq);
           dq.pop_front();
           }
           n++;
           }
           */
    }

};

Array<Model, 16*4> &bires_models() {
    static BiResModels ms;
    return ms.models;
}



END_JN

