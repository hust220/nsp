#pragma once

#include "pdb.hpp"
#include "geom.hpp"
#include "env.hpp"

namespace jian {

class BuildLine {
public:
    Chain m_chain;
    Chain m_frag;
    Residue m_res;

    BuildLine() {
        load_res();
    }

    void load_res() {
        chain_read_model(m_frag, to_str(Env::lib(),"/RNA/pars/nuc3d/BuildChain/frag.pdb"));
        m_res = m_frag[0];
    }

    BuildLine &operator ()(int n) {
        m_chain.clear();
        for (int i = 0; i < n; i++) {
            m_chain.push_back(m_res);
            for (auto && at : m_res) {
                at[0] += 10;
            }
        }
        return *this;
    }
};


}


