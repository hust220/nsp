#include "rna_mc.hpp"
#include "rna_mc_backup.hpp"
#include "rna_mc_select.hpp"

namespace jian {

void dhmc_backup(DHMC &m) {
    m._moved_atoms.clear();

    Int l = size(m._seq);

    for (int i = 0; i < l; i++) {
        if (m.is_selected(i)) {
            for (auto && atom : m._pred_chain[i]) {
                m._moved_atoms.push_back(atom);
            }
        }
    }
}

}
