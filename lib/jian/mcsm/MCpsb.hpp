#pragma once

#include "MCxp.hpp"
#include "MCen.hpp"
#include "../cg.hpp"
#include "../scoring/ScorePsb.hpp"

#define MEM_EN_MCPSB len, ang, dih, crash, cons, vdw, stacking, pairing
#define DEF_MEM_EN_MCPSB(a) double a = 0;
#define SUM_MEM_EN_MCPSB(a) + a
#define PRINT_MEM_EN_MCPSB(a) << a << PP_STRING3((a)) << ' '
#define PRINT_MEM_MCPSB(a) << e.a << PP_STRING3((a)) << ' '

#define MCPSB_ENERGY_CRASH(name, cond1, cond2) \
    void mc_##name##_energy_crash(en_t &e) { \
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
									e.crash += _mc_crash_weight * en_crash(_pred_chain[p.first], _pred_chain[p.second]); \
									e.pairing += _mc_pairing_weight * m_scorer.en_pairing(_pred_chain[p.first], _pred_chain[p.second]);\
								} \
							} \
						} \
					} \
				} \
			} \
		} \
	}

#define MCPSB_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d; \
        for (auto && n : m_continuous_pts) { \
             cond { \
                d = geom::distance(_pred_chain[n][0], _pred_chain[n+1][0]); \
                e.len += _mc_bond_length_weight * square(d - 6.1); \
				e.crash += _mc_crash_weight * en_crash(_pred_chain[n], _pred_chain[n+1]);\
				e.pairing += _mc_pairing_weight * m_scorer.en_pairing(_pred_chain[n], _pred_chain[n+1]);\
             } \
        } \
    }

#define MCPSB_ENERGY_ANGLE(name, cond) \
    void mc_##name##_energy_angle(en_t &e) { \
        double d; \
        int len = _seq.size(); \
        for (auto && i : m_ang_pts) { \
             cond { \
                d = geom::angle(_pred_chain[i][0],  \
                                _pred_chain[i+1][0],  \
                                _pred_chain[i+2][0]); \
                e.ang += _mc_bond_angle_weight * square(d - _mc_bond_angle_std); \
            } \
        } \
    }

#define MCPSB_ENERGY_DIHEDRAL(name, cond) \
    void mc_##name##_energy_dihedral(en_t &e) { \
        double d; \
        int len = _seq.size(); \
        for (auto && i : m_dih_pts) { \
            cond { \
                d = geom::dihedral(_pred_chain[i][0],  \
                                   _pred_chain[i+1][0],  \
                                   _pred_chain[i+2][0],  \
                                   _pred_chain[i+3][0]); \
                d = d - _mc_bond_dihedral_std; \
                d = 3.3 - 4 * std::cos(d) + std::cos(2 * d) - 0.44 * std::cos(3 * d); \
                e.dih += _mc_bond_dihedral_weight * d; \
            } \
        } \
    }


#define MCPSB_ENERGY_CONSTRAINTS(name, cond) \
    void mc_##name##_energy_constraints(en_t &e) { \
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

			template<>
			class MCen<CGpsb> : public MCxp<CGpsb> {
			public:
				struct en_t {
					JN_MAP(DEF_MEM_EN_MCPSB, MEM_EN_MCPSB);
					double sum() const { return 0 JN_MAP(SUM_MEM_EN_MCPSB, MEM_EN_MCPSB); }
					void print() const { LOG << sum() << "(total) " JN_MAP(PRINT_MEM_EN_MCPSB, MEM_EN_MCPSB) << std::endl; }
				};

				//Mat m_stacking_pars;
				//Mat m_pairing_pars;
				std::vector<int> m_indices;
				Score<CGpsb> m_scorer;

				MCen() = default;

				void init(const Par &par) {
					MCxp<CGpsb>::init(par);

					LOG << "# Initializing scorer..." << std::endl;
					m_scorer.init();

					//LOG << "# Reading stacking and pairing parameters..." << std::endl;
					//read_stacking_pairing_parameters();

					LOG << "# Seting indices..." << std::endl;
					set_indices();
				}

				void set_indices() {
					std::map<char, int> m{ {'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3} };
					int i;
					m_indices.resize(_seq.size());
					for (i = 0; i < _seq.size(); i++) {
						m_indices[i] = m[_seq[i]];
					}
				}

				void print_pairing() {
					Mat arr(3, 3);
					double d;
					int i, j;

					for (i = 0; i < _seq.size(); i++) {
						for (j = i + 1; j < _seq.size(); j++) {
							//set_arr(arr, i, j);
							d = m_scorer.en_pairing(_pred_chain[i], _pred_chain[j]);
							if (d != 0) {
								LOG << i + 1 << ' ' << j + 1 << ' ' << d << std::endl;
							}
						}
					}
				}

				virtual double mc_partial_energy() {
					en_t e;
					mc_partial_energy_crash(e);
					mc_partial_energy_bond(e);
					mc_partial_energy_angle(e);
					mc_partial_energy_dihedral(e);
					mc_partial_energy_constraints(e);
					return e.sum();
				}

				void mc_total_energy(en_t &e) {
					mc_total_energy_crash(e);
					mc_total_energy_bond(e);
					mc_total_energy_angle(e);
					mc_total_energy_dihedral(e);
					mc_total_energy_constraints(e);
				}

				MCPSB_ENERGY_CRASH(partial, if (is_selected(n)), if (!(is_selected(t)) && (t - n != 1 && n - t != 1)));
				MCPSB_ENERGY_CRASH(total, , if (t - n > 1));

				MCPSB_ENERGY_BOND(partial, if (is_selected(n) ^ is_selected(n + 1)));
				MCPSB_ENERGY_BOND(total, );

				MCPSB_ENERGY_ANGLE(partial, if ((is_selected(i) + is_selected(i + 1) + is_selected(i + 2)) % 3 != 0));
				MCPSB_ENERGY_ANGLE(total, );

				MCPSB_ENERGY_DIHEDRAL(partial, if ((is_selected(i) + is_selected(i + 1) + is_selected(i + 2) + is_selected(i + 3)) % 4 != 0));
				MCPSB_ENERGY_DIHEDRAL(total, );

				MCPSB_ENERGY_CONSTRAINTS(partial, if (is_selected(c.key[0]) ^ is_selected(c.key[1])));
				MCPSB_ENERGY_CONSTRAINTS(total, );

				double en_crash(const Residue &r1, const Residue &r2) {
					int i, j;
					double d, e;

					e = 0;
					for (i = 0; i < 3; i++) {
						for (j = 0; j < 3; j++) {
							d = geom::distance(r1[i], r2[j]);
							if (i == 0 && j == 0) {
								if (d < 5) {
									e += square(d - 5);
								}
							}
							else if (i == 1 && j == 1) {
								if (d < 5) {
									e += square(d - 5);
								}
							}
							else if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
								if (d < 6.5) {
									e += square(d - 6.5);
								}
							}
							else if (d < 3.5) {
								e += square(d - 3.5);
							}
						}
					}
					return e;
				}

				virtual double dist_two_res(const Residue &r1, const Residue &r2) const {
					return geom::distance(r1[0], r2[0]);
				}

				double total_energy() {
					en_t e;
					mc_total_energy(e);
					e.print();
					return 0;
				}

				virtual void write_en() {
					en_t e;
					mc_total_energy(e);
					LOG << _mc_step + 1 << ": " << e.sum() << "(total) "
						JN_MAP(PRINT_MEM_MCPSB, MEM_EN_MCPSB) << _mc_tempr << "(tempr) "
						<< _mc_local_succ_rate << "(rate)" << std::endl;
				}

				virtual std::string file_parameters() const {
					return "mcpsb";
				}

				virtual void finish_run() {
					//LOG << "# Displaying stacking information" << std::endl;
					//print_stacking();

					LOG << "# Displaying pairing information" << std::endl;
					print_pairing();
				}

				virtual void mc_select() = 0;
				virtual bool is_selected(const int &i) const = 0;
				virtual Vec rotating_center() const = 0;

			};

			using MCpsb = MCen<CGpsb>;

		} // namespace mc
	} // namespace nuc3d
} // namespace jian

