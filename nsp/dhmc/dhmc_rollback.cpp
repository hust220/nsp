#include "DHMC.hpp"
#include "dhmc_rollback.hpp"
#include "dhmc_select.hpp"

BEGIN_JN

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

END_JN
