#pragma once

#include "pdb.hpp"
#include "geom.hpp"
#include "env.hpp"

BEGIN_JN

class BuildChain {
public:
    Chain m_chain;
    Chain m_frag;
    Mat m_mat_chain;
    Mat m_mat_frag;
    int m_len;
    std::vector<Str> m_joints;

    BuildChain() : m_joints{"O5*", "C3*", "C1*"} {
        load_frag();
        set_mat(m_mat_frag, m_frag, 0);
    }

    void load_frag() {
        chain_read_model(m_frag, to_str(Env::lib(),"/RNA/pars/nuc3d/BuildChain/frag.pdb"));
        m_len = size(m_frag);
    }

    void set_mat(Mat &mat, const Chain &chain, int beg) {
        int n;
        mat.resize(3*(m_len-1), 3);
        for (int i = 0; i < m_len-1; i++) {
            for (auto && atom : chain[i+beg]) {
                auto it = std::find(m_joints.begin(), m_joints.end(), atom.name);
                if (it != m_joints.end()) {
                    n = std::distance(m_joints.begin(), it);
                    for (int j = 0; j < 3; j++) {
                        mat(i * 3 + n, j) = atom[j];
                    }
                }
            }
        }
    }

    void increment() {
        set_mat(m_mat_chain, m_chain, size(m_chain) - m_len + 1);
        Residue r = m_frag[m_len-1];
        geom::Superposition<Num> sp(m_mat_frag, m_mat_chain);
        for (auto && atom : r) {
            sp.apply(atom);
        }
        m_chain.push_back(std::move(r));
    }

    BuildChain &operator ()(int n) {
        if (n < m_len) {
            throw std::invalid_argument(to_str("Can't build a chain less than ", m_len, " nt!"));
        }
        m_chain = m_frag;
        for (int i = 0; i < n - m_len; i++) {
            increment();
        }
        return *this;
    }
};


}


