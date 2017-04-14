#include "DHMC.hpp"
#include "dhmc_backup.hpp"
#include "dhmc_select.hpp"

BEGIN_JN

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

END_JN
