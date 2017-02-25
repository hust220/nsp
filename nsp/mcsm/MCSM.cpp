#include "MCSM.hpp"

#define PRINT_MEM_MCPSB(a) << e.a << PP_STRING3((a)) << ' '

BEGIN_JN

void MCSM::init(const Par &par) {
    MCBase::init(par);

    log << "# Initializing scorer..." << std::endl;
    m_scorer = Score::fac_t::create(m_cg_type);
    m_scorer->init();

    log << "# Seting indices..." << std::endl;
    set_indices();
}

void MCSM::set_indices() {
    std::map<char, int> m{ { 'A', 0 },{ 'U', 1 },{ 'T', 1 },{ 'G', 2 },{ 'C', 3 } };
    int i;
    m_indices.resize(_seq.size());
    for (i = 0; i < _seq.size(); i++) {
        m_indices[i] = m[_seq[i]];
    }
}

void MCSM::print_pairing() {
    Mat arr(3, 3);
    int i, j;

    for (i = 0; i < _seq.size(); i++) {
        for (j = i + 1; j < _seq.size(); j++) {
            m_scorer->en_bp(_pred_chain[i], _pred_chain[j]);
            if (m_scorer->m_en_pairing != 0 || m_scorer->m_en_stacking != 0) {
                log
                    << i + 1 << ' ' << j + 1 << ' '
                    << m_scorer->m_en_stacking << ' '
                    << m_scorer->m_en_pairing << ' '
                    //<< m_scorer->m_en_wc << ' '
                    //<< m_scorer->m_en_nwc << ' '
                    << std::endl;
            }
        }
    }
}

double MCSM::mc_partial_energy() {
    en_t e;
    mc_energy_crash(e, false);
    mc_energy_bond(e, false);
    mc_energy_angle(e, false);
    mc_energy_dihedral(e, false);
    mc_energy_constraints(e, false);
    return e.sum();
}

void MCSM::set_total_energy(en_t &e) {
    mc_energy_crash(e, true);
    mc_energy_bond(e, true);
    mc_energy_angle(e, true);
    mc_energy_dihedral(e, true);
    mc_energy_constraints(e, true);
}

double MCSM::mc_total_energy() {
    en_t e;
    set_total_energy(e);
    return e.sum();
}

void MCSM::mc_energy_crash(en_t &e, bool is_total) {
    int a, b, c, i, j, k, n;
    Num d;
    for (n = 0; n < _seq.size(); n++) {
        if (is_total || is_selected(n)) {
            for (i = -m_box; i <= m_box; i++) for (j = -m_box; j <= m_box; j++) for (k = -m_box; k <= m_box; k++) {
                SpaceItem &it = item(n);
                a = space_index(it[0]) + i;
                b = space_index(it[1]) + j;
                c = space_index(it[2]) + k;
                SpaceVal &s = m_space[a][b][c];
                for (auto && t : s) {
                    if ((is_total && t - n > 0) || (!is_total && ((is_selected(t) && t - n > 0) || (!is_selected(t))))) {
                    //if ((is_total && t - n > 0) || (!is_total && !is_selected(t))) {
                        auto p = std::minmax(n, t);
                        e.crash += _mc_crash_weight * m_scorer->en_crash(_pred_chain[p.first], _pred_chain[p.second]);
                        if (!m_grow_mode || p.first < m_grow_length && p.second < m_grow_length) {
                            m_scorer->en_bp(_pred_chain[p.first], _pred_chain[p.second]);
                            d = _mc_pairing_weight * m_scorer->m_en_pairing;
                            if (m_bps[p.first] == p.second) e.cons += _mc_constraints_weight * d;
                            e.pairing += d/* * (m_bps[p.first] == p.second ? 10 : 1)*/;
                            e.stacking += _mc_stacking_weight * m_scorer->m_en_stacking;
                            e.vdw += _mc_vdw_weight * m_scorer->m_en_vdw;
                        }
                    }
                }
            }
        }
    }
}

void MCSM::mc_energy_bond(en_t &e, bool is_total) {
    for (auto && n : m_continuous_pts) {
        if (!m_grow_mode || n < m_grow_length) {
            if (is_total || is_selected(n) || is_selected(n + 1)) {
            //if (is_total || (is_selected(n) + is_selected(n + 1)) % 2 != 0) {
                e.len += _mc_bond_length_weight * m_scorer->en_len(_pred_chain, n);
            }
        }
    }
}

void MCSM::mc_energy_angle(en_t &e, bool is_total) {
    int len = _seq.size();
    for (auto && i : m_ang_pts) {
        if (!m_grow_mode || i < m_grow_length) {
            if (is_total || is_selected(i) || is_selected(i + 1) || is_selected(i + 2)) {
            //if (is_total || (is_selected(i) + is_selected(i + 1) + is_selected(i + 2)) % 3 != 0) {
                e.ang += _mc_bond_angle_weight * m_scorer->en_ang(_pred_chain, i);
            }
        }
    }
}

void MCSM::mc_energy_dihedral(en_t &e, bool is_total) {
    int len = _seq.size();
    for (auto && i : m_dih_pts) {
        if (!m_grow_mode || i < m_grow_length) {
            if (is_total || is_selected(i) || is_selected(i + 1) || is_selected(i + 2) || is_selected(i + 3)) {
            //if (is_total || (is_selected(i) + is_selected(i + 1) + is_selected(i + 2) + is_selected(i + 3)) % 4 != 0) {
                e.dih += _mc_bond_dihedral_weight * m_scorer->en_dih(_pred_chain, i);
            }
        }
    }
}

void MCSM::mc_energy_constraints(en_t &e, bool is_total) {
    double d, k;

    if (m_cal_en_constraints) {
        // residue contacts
        k = _mc_constraints_weight / _constraints.contacts.size();
        for (auto && c : _constraints.contacts) {
            if (!m_grow_mode || c.key[0] < m_grow_length && c.key[1] < m_grow_length) {
                if (is_total || is_selected(c.key[0]) || is_selected(c.key[1])) {
                //if (is_total || is_selected(c.key[0]) ^ is_selected(c.key[1])) {
                    d = geom::distance(_pred_chain[c.key[0]][0], _pred_chain[c.key[1]][0]);
                    if (d < 17) {
                        e.cons += -100.0 * k;
                    }
                    else {
                        e.cons += square(d - 17) * k;
                    }
                }
            }
        }

        // residue distances
        k = _mc_constraints_weight / _constraints.distances.size();
        for (auto && c : _constraints.distances) {
            if (!m_grow_mode || c.key[0] < m_grow_length && c.key[1] < m_grow_length) {
                if (is_total || is_selected(c.key[0]) || is_selected(c.key[1])) {
                //if (is_total || is_selected(c.key[0]) ^ is_selected(c.key[1])) {
                    d = geom::distance(_pred_chain[c.key[0]][0], _pred_chain[c.key[1]][0]);
                    e.cons += square(d - c.value) * k;
                }
            }
        }

        // base pairing
        k = _mc_constraints_weight / m_distance_constraints.size();
        for (auto && c : m_distance_constraints) {
            if (!m_grow_mode || c.atom1.n_res < m_grow_length && c.atom2.n_res < m_grow_length) {
                if (is_total || is_selected(c.atom1.n_res) || is_selected(c.atom2.n_res)) {
                //if (is_total || is_selected(c.atom1.n_res) ^ is_selected(c.atom2.n_res)) {
                    d = geom::distance(_pred_chain[c.atom1.n_res][c.atom1.n_atom], _pred_chain[c.atom2.n_res][c.atom2.n_atom]);
                    e.cons += (d < c.min ? square(d - c.min) : (d > c.max ? square(d - c.max) : 0)) * k;
                }
            }
        }
    }
}

double MCSM::dist_two_res(const Residue &r1, const Residue &r2) const {
    return geom::distance(r1[0], r2[0]);
}

void MCSM::write_en() {
    en_t e;
    set_total_energy(e);
    log 
        << _mc_step + 1 << ": " << e.sum() << "(total) "
        JN_MAP(PRINT_MEM_MCPSB, MEM_EN_MCPSB)
        << _mc_tempr << "(tempr) "
        << _mc_local_succ_rate << "(rate)" << std::endl;
}

S MCSM::file_parameters() const {
    return "mcpsb";
}

void MCSM::finish_run() {
    log << "# Displaying pairing information" << std::endl;
    print_pairing();
}

END_JN
