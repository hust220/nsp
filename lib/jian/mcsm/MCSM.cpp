#include "MCSM.hpp"

#define PRINT_MEM_MCPSB(a) << e.a << PP_STRING3((a)) << ' '

#define MCPSB_ENERGY_CRASH(name, cond1, cond2) \
    void MCSM::mc_##name##_energy_crash(en_t &e) { \
        int a, b, c, i, j, k, n; \
        for (n = 0; n < _seq.size(); n++) { \
            cond1 { \
                for (i = -m_box; i <= m_box; i++) { \
                    for (j = -m_box; j <= m_box; j++) { \
                        for (k = -m_box; k <= m_box; k++) { \
                            item_t &it = item(n); \
                            a = space_index(it[0])+i; \
                            b = space_index(it[1])+j; \
                            c = space_index(it[2])+k; \
                            space_val_t &s = m_space[a][b][c]; \
                            for (auto && t : s) { \
                                cond2 { \
                                    auto p = std::minmax(n, t); \
									e.crash += _mc_crash_weight * m_scorer->en_crash(_pred_chain[p.first], _pred_chain[p.second]); \
									m_scorer->en_bp(_pred_chain[p.first], _pred_chain[p.second]);\
									/*e.pairing += _mc_pairing_weight * m_scorer->m_en_pairing;*/\
									e.wc += _mc_wc_weight * m_scorer->m_en_wc;\
									e.nwc += _mc_nwc_weight * m_scorer->m_en_nwc;\
									e.stacking += _mc_stacking_weight * m_scorer->m_en_stacking;\
								} \
							} \
						} \
					} \
				} \
			} \
		} \
	}

#define MCPSB_ENERGY_BOND(name, cond) \
    void MCSM::mc_##name##_energy_bond(en_t &e) { \
        for (auto && n : m_continuous_pts) { \
             cond { \
				e.len += _mc_bond_length_weight * m_scorer->en_len(_pred_chain[n], _pred_chain[n+1]);\
				e.crash += _mc_crash_weight * m_scorer->en_crash(_pred_chain[n], _pred_chain[n+1]);\
				m_scorer->en_bp(_pred_chain[n], _pred_chain[n+1]);\
				e.stacking += _mc_stacking_weight * m_scorer->m_en_stacking;\
             } \
        } \
    }

#define MCPSB_ENERGY_ANGLE(name, cond) \
    void MCSM::mc_##name##_energy_angle(en_t &e) { \
        int len = _seq.size(); \
        for (auto && i : m_ang_pts) { \
             cond { \
				e.ang += _mc_bond_angle_weight * m_scorer->en_ang(_pred_chain[i], _pred_chain[i+1], _pred_chain[i+2]);\
            } \
        } \
    }

#define MCPSB_ENERGY_DIHEDRAL(name, cond) \
    void MCSM::mc_##name##_energy_dihedral(en_t &e) { \
        int len = _seq.size(); \
        for (auto && i : m_dih_pts) { \
            cond { \
				e.dih += _mc_bond_angle_weight * m_scorer->en_dih(_pred_chain[i], _pred_chain[i+1], _pred_chain[i+2], _pred_chain[i+3]);\
            } \
        } \
    }


#define MCPSB_ENERGY_CONSTRAINTS(name, cond) \
    void MCSM::mc_##name##_energy_constraints(en_t &e) { \
        double d, k; \
		k = _mc_constraints_weight / _constraints.contacts.size(); \
		for (auto && c : _constraints.contacts) { \
			cond { \
				d = geom::distance(_pred_chain[c.key[0]][0], _pred_chain[c.key[1]][0]); \
				if (d < 17) { \
					e.cons += -100.0 * k; \
				} else { \
					e.cons += square(d - 17) * k; \
				} \
			} \
		}\
		k = _mc_constraints_weight / _constraints.distances.size(); \
		for (auto && c : _constraints.distances) { \
			cond { \
				d = geom::distance(_pred_chain[c.key[0]][0], _pred_chain[c.key[1]][0]); \
				e.cons += square(d - c.value) * k; \
			} \
		}\
	}

namespace jian {
	namespace nuc3d {
		namespace mc {

			void MCSM::init(const Par &par) {
				MCBase::init(par);

				LOG << "# Initializing scorer..." << std::endl;
				m_scorer = ScoreBase::fac_t::create(m_cg_type);
				m_scorer->init();

				LOG << "# Seting indices..." << std::endl;
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
							LOG 
								<< i + 1 << ' ' << j + 1 << ' ' 
								<< m_scorer->m_en_stacking << ' ' 
								<< m_scorer->m_en_wc << ' ' 
								<< m_scorer->m_en_nwc << std::endl;
						}
					}
				}
			}

			double MCSM::mc_partial_energy() {
				en_t e;
				mc_partial_energy_crash(e);
				mc_partial_energy_bond(e);
				mc_partial_energy_angle(e);
				mc_partial_energy_dihedral(e);
				mc_partial_energy_constraints(e);
				return e.sum();
			}

			void MCSM::mc_total_energy(en_t &e) {
				mc_total_energy_crash(e);
				mc_total_energy_bond(e);
				mc_total_energy_angle(e);
				mc_total_energy_dihedral(e);
				mc_total_energy_constraints(e);
			}

			MCPSB_ENERGY_CRASH(partial, if (is_selected(n)), if (!(is_selected(t)) && (t - n != 1 && n - t != 1)));
			MCPSB_ENERGY_CRASH(total, , if (t - n > 1));

			MCPSB_ENERGY_BOND(partial, if (is_selected(n) || is_selected(n + 1)));
			MCPSB_ENERGY_BOND(total, );

			MCPSB_ENERGY_ANGLE(partial, if ((is_selected(i) + is_selected(i + 1) + is_selected(i + 2)) % 3 != 0));
			MCPSB_ENERGY_ANGLE(total, );

			MCPSB_ENERGY_DIHEDRAL(partial, if ((is_selected(i) + is_selected(i + 1) + is_selected(i + 2) + is_selected(i + 3)) % 4 != 0));
			MCPSB_ENERGY_DIHEDRAL(total, );

			MCPSB_ENERGY_CONSTRAINTS(partial, if (is_selected(c.key[0]) ^ is_selected(c.key[1])));
			MCPSB_ENERGY_CONSTRAINTS(total, );

			double MCSM::dist_two_res(const Residue &r1, const Residue &r2) const {
				return geom::distance(r1[0], r2[0]);
			}

			double MCSM::total_energy() {
				en_t e;
				mc_total_energy(e);
				e.print();
				return 0;
			}

			void MCSM::write_en() {
				en_t e;
				mc_total_energy(e);
				LOG << _mc_step + 1 << ": " << e.sum() << "(total) "
					JN_MAP(PRINT_MEM_MCPSB, MEM_EN_MCPSB)
					<< _mc_tempr << "(tempr) "
					<< _mc_local_succ_rate << "(rate)" << std::endl;
			}

			std::string MCSM::file_parameters() const {
				return "mcpsb";
			}

			void MCSM::finish_run() {
				LOG << "# Displaying pairing information" << std::endl;
				print_pairing();
			}

		}
	}
} // namespace jian
