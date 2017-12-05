#include "nsp.hpp"
#include "score_par_bp.hpp"
#include "dca.hpp"
#include "rss_get_ss.hpp"

BEGIN_JN

namespace {
    REGISTER_NSP_COMPONENT(ss) {
        auto g = par.getv("global");
        Str molfile = g[1];
        Model m;
        mol_read(m, molfile);
        JN_OUT << get_ss(m.residues()) << std::endl;
        /*
           Par::pars_t mols;

           par.setv(mols, "mols");
           for (auto && mol : mols) {
           for_each_model(mol, [&mol](const Model &model, int i) {
           JN_OUT << to_str(mol, ":model-", i + 1) << std::endl;
           JN_OUT << get_ss(model) << std::endl;
           });
           }
           */
    }
}
END_JN
















