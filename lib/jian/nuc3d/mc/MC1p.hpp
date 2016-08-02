#pragma once

#include "MCxp.hpp"
#include "../../cg.hpp"

namespace jian {
namespace nuc3d {
namespace mc {

class MC1p : public MCxp<CG1p> {
public:    
    #define MC1p_en_t_m len, ang, dih, crash, cons, vdw
    struct en_t {
        #define MC1p_en_t_m_def(a) double a = 0;
        JN_MAP(MC1p_en_t_m_def, MC1p_en_t_m)
        #define MC1p_en_t_m_sum(a) + a
        double sum() const { return 0 JN_MAP(MC1p_en_t_m_sum, MC1p_en_t_m); }
    };

    MC1p(const Par &par) : MCxp(par) {}

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

#define MC1P_ENERGY_CRASH(name, cond1, cond2) \
    void mc_##name##_energy_crash(en_t &e) { \
        int a, b, c; \
        double d; \
        for (int n = 0; n < _seq.size(); n++) { \
            cond1 { \
                for (int i = -m_box; i <= m_box; i++) { \
                    for (int j = -m_box; j <= m_box; j++) { \
                        for (int k = -m_box; k <= m_box; k++) { \
                            item_t &it = item(n); \
                            a = space_index(it[0])+i; \
                            b = space_index(it[1])+j; \
                            c = space_index(it[2])+k; \
                            space_val_t &s = m_space[a][b][c]; \
                            for (auto && t : s) { \
                                cond2 { \
                                    d = geom::distance(_pred_chain[n][0], _pred_chain[t][0]); \
                                    if (d < 9) { \
                                        e.crash += _mc_crash_weight * square(d - 9); \
                                    } \
                                } \
                            } \
                        } \
                    } \
                } \
            } \
        } \
    }

    MC1P_ENERGY_CRASH(partial, if (is_selected(n)), if (!(is_selected(t)) && (t - n != 1 && n - t != 1)))
    MC1P_ENERGY_CRASH(total, , if (t - n > 1))

#define MC1P_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d; \
        Mat arr = Mat::Zero(3, 3); \
        for (auto && n : m_continuous_pts) { \
             cond { \
                d = geom::distance(_pred_chain[n][0], _pred_chain[n+1][0]); \
                e.len += _mc_bond_length_weight * square(d - 6.1); \
             } \
        } \
    }
    MC1P_ENERGY_BOND(partial, if (is_selected(n) ^ is_selected(n+1)))
    MC1P_ENERGY_BOND(total, )

#define MC1P_ENERGY_ANGLE(name, cond) \
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
    MC1P_ENERGY_ANGLE(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2)) % 3 != 0))
    MC1P_ENERGY_ANGLE(total, )

#define MC1P_ENERGY_DIHEDRAL(name, cond) \
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
    MC1P_ENERGY_DIHEDRAL(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2) + is_selected(i+3)) % 4 != 0))
    MC1P_ENERGY_DIHEDRAL(total, )

#define MC1P_ENERGY_CONSTRAINTS(name, cond) \
    void mc_##name##_energy_constraints(en_t &e) { \
        double d; \
        int n = _constraints.distances.size(); \
        for (auto && row : _constraints.distances) { \
            cond { \
                d = geom::distance(_pred_chain[row.key[0]][0], _pred_chain[row.key[1]][0]); \
                if (d < 17) { \
                    e.cons += _mc_constraints_weight * (-100); \
                } else { \
                    e.cons += _mc_constraints_weight * 0.1 * row.weight * square(d - 17); \
                } \
            } \
        } \
    }
    MC1P_ENERGY_CONSTRAINTS(partial,if (is_selected(row.key[0]) ^ is_selected(row.key[1])))
    MC1P_ENERGY_CONSTRAINTS(total,)

    virtual double dist_two_res(const Residue &r1, const Residue &r2) const {
        return geom::distance(r1[0], r2[0]);
    }

    virtual void write_en() {
        en_t e;
        mc_total_energy(e);
        #define MC1p_print(a) << e.a << PP_STRING3((a)) << ' '
        LOG << _mc_step + 1 << ": " <<  e.sum() << "(total) "
               JN_MAP(MC1p_print, MC1p_en_t_m) << _mc_tempr << "(tempr) " << _mc_local_succ_rate << "(rate)"  << std::endl;
    }

    virtual std::string file_parameters() const {
        return "mc1p";
    }

    virtual void mc_select() = 0;
    virtual bool is_selected(const int &i) const = 0;
    virtual Vec rotating_center() const = 0;

};

} // namespace mc
} // namespace nuc3d
} // namespace jian

