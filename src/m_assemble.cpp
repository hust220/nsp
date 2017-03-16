#include "nsp.hpp"
#include <nsp/nuc3d/Assemble.hpp>
#include <nsp/pdb/utils/cluster_chains.hpp>
#include <nsp/scoring/Score.hpp>

BEGIN_JN

namespace {

    void write_pred(const nuc3d::Assemble &ass, int n) {
        mol_write(ass._pred_chain, to_str(ass.m_out_dir, '/', ass._name, ".pred", n, ".pdb"));
    }

    REGISTER_NSP_COMPONENT(assemble) {
        nuc3d::Assemble ass(par);

        std::ostringstream stream;
        int n;

        int num = 1;
        par.set(num, "n", "num");

        ass.select_templates();
        ass.assemble();

        n = 1;
        write_pred(ass, n);
        //mol_write(ass._pred_chain, to_str(ass._name, ".", n, ".pred.pdb"));
        for (n = 2; n <= num; n++) {
            ass.sample_all_templates();
            ass.assemble();
            ass.log << "# Writing sampling structure " << n << std::endl;
            write_pred(ass, n);
            //mol_write(ass._pred_chain, to_str(ass._name, ".pred", n, ".pdb"));
        }

    }

}

END_JN

