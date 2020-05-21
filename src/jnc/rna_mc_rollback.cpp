#include "rna_mc.hpp"
#include "rna_mc_rollback.hpp"
#include "rna_mc_select.hpp"

namespace jian {

void dhmc_rollback(DHMC &m) {
    Int l = size(m._seq);

    for (int i = 0; i < l; i++) {
        if (m.is_selected(i)) {
            for (auto && atom : m._pred_chain[i]) {
                atom = m._moved_atoms.front();
                m._moved_atoms.pop_front();
            }
            m.space_update_item(i);
        }
    }
}

}
