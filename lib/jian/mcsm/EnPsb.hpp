#pragma once

#include "MCxp.hpp"
#include "MCen.hpp"
#include "../cg.hpp"

namespace jian {
namespace nuc3d {
namespace mc {

template<>
class En<CGpsb> {
public:    
    #define EnPsb_en_t_m len, ang, dih, crash, cons, vdw, stacking, pairing
    struct en_t {
        #define EnPsb_en_t_m_def(a) double a = 0;
        JN_MAP(EnPsb_en_t_m_def, EnPsb_en_t_m)
        #define EnPsb_en_t_m_sum(a) + a
        double sum() const { return 0 JN_MAP(EnPsb_en_t_m_sum, EnPsb_en_t_m); }
        #define EnPsb_print_en(a) << a << PP_STRING3((a)) << ' '
        void print() const {LOG << sum() << "(total) " JN_MAP(EnPsb_print_en, EnPsb_en_t_m) << std::endl;}
        #undef EnPsb_print_en
    };

    Mat m_stacking_pars;
    Mat m_pairing_pars;

    void init() {
        read_stacking_pairing_parameters();
    }

    void read_stacking_pairing_parameters() {
        m_stacking_pars.resize(8, 8);
        std::string file_name = Env::lib() + "/RNA/pars/nuc3d/mc/mcpsb.stacking.par";
        std::ifstream ifile(file_name.c_str());
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) ifile >> m_stacking_pars(i, j);
        ifile.close();
        m_pairing_pars.resize(3, 3);
        file_name = Env::lib() + "/RNA/pars/nuc3d/mc/mcpsb.pairing.par";
        ifile.open(file_name.c_str());
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) ifile >> m_pairing_pars(i, j);
        ifile.close();
    }

    void set_arr(Mat &arr, int m, int n) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                arr(i, j) = geom::distance(_pred_chain[m][i], _pred_chain[n][j]);
            }
        }
    }

    void print_stacking() {
        Mat arr(3, 3);
        double d;
        for (int i = 0; i < _seq.size(); i++) {
            for (int j = i + 1; j < _seq.size(); j++) {
                set_arr(arr, i, j);
                d = en_stacking(arr, i, j);
                if (d != 0) {
                    LOG << i+1 << ' ' << j+1 << ' ' << d  << std::endl;
                }
            }
        }
    }

    void print_pairing() {
        Mat arr(3, 3);
        double d;
        for (int i = 0; i < _seq.size(); i++) {
            for (int j = i + 1; j < _seq.size(); j++) {
                set_arr(arr, i, j);
                d = en_pairing(arr, i, j);
                if (d != 0) {
                    LOG << i+1 << ' ' << j+1 << ' ' << d  << std::endl;
                }
            }
        }
    }

    double en_stacking(const Mat &m, int a, int b) const {
        double en = 0, d;
        int x = m_indices[a] * 2, y = m_indices[b] * 2;
        if (m(0, 0) > 4 && m(0, 0) < 7 && m(1, 1) > 3.5 && m(1, 1) < 5.5) {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    d = square(m(i, j) - m_stacking_pars(x + i, y + j));
                    if (d > 1) return 0;
                    en += d;
                }
            }
            return (-40 + en) * _mc_stacking_weight;
        } else {
            return 0;
        }
    }

    double en_pairing(const Mat & m, int a, int b) const {
        double d, en;
        int x, y;

        en = 0;
        x = a+b;
        y = a*b;
        if (m(0, 0) > 14.5 && m(0, 0) < 17.5 && m(1, 1) > 9.5 && m(1, 1) < 14.5) {
            if (a == 0 || a == 2) {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        d = square(m(i, j) - m_pairing_pars(i, j));
                        if (d > 1) return 0;
                        en += d;
                    }
                }
            } else {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        d = square(m(i, j) - m_pairing_pars(j, i));
                        if (d > 1) return 0;
                        en += d;
                    }
                }
            }
            if (x == 1) {
                return (-100+en) * 2 * m_pairing_weight;
            } else if (x == 5) {
                return (-100+en) * 3 * m_pairing_weight;
            } else if (y == 2) {
                return (-100+en) * 1 * m_pairing_weight;
            } else {
                return (-100+en) * 0.5 * m_pairing_weight;
            }
        } else {
            return 0;
        }
    }

    double en_crash(const Residue & a, const Residue & b, Mat & arr) {
        double d, e;
        unsigned int i, j;

        e = 0;
        for (i = 0; i < 3; i++) { 
            for (j = 0; j < 3; j++) { 
                d = geom::distance(a[i], b[j]); 
                arr(i, j) = d; 
                if (i == 0 || j == 0) { 
                    if (d < 7) { 
                        e += square(d - 7); 
                    } 
                } else if ((i == 1 || j == 1) && d < 5) { 
                    e += square(d - 5); 
                } else if (d < 3.5) { 
                    e += square(d - 3.5); 
                } 
            } 
        } 
        return e * m_crash_weight;
    }

#define MCPSB_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d; \
        Mat arr = Mat::Zero(3, 3); \
        for (auto && n : m_continuous_pts) { \
             cond { \
                d = geom::distance(_pred_chain[n][0], _pred_chain[n+1][0]); \
                e.len += _mc_bond_length_weight * square(d - 6.1); \
                for (int i = 0; i < 3; i++) {  \
                    for (int j = 0; j < 3; j++) {  \
                        d = geom::distance(_pred_chain[n][i], _pred_chain[n+1][j]);  \
                        arr(i, j) = d;  \
                        if (d < 3.5) {  \
                            e.crash += _mc_crash_weight * square(d - 3.5);  \
                        }  \
                    }  \
                }  \
                e.stacking += en_stacking(arr, n, n+1);  \
             } \
        } \
    }
    MCPSB_ENERGY_BOND(partial, if (is_selected(n) ^ is_selected(n+1)))
    MCPSB_ENERGY_BOND(total, )

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
    MCPSB_ENERGY_ANGLE(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2)) % 3 != 0))
    MCPSB_ENERGY_ANGLE(total, )

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
    MCPSB_ENERGY_DIHEDRAL(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2) + is_selected(i+3)) % 4 != 0))
    MCPSB_ENERGY_DIHEDRAL(total, )

#define MCPSB_ENERGY_CONSTRAINTS(name, cond) \
    void mc_##name##_energy_constraints(en_t &e) { \
        double d; \
        int n = _constraints.distances.size(); \
        for (auto && row : _constraints.distances) { \
            cond { \
                d = geom::distance(_pred_chain[row.key[0]][0], _pred_chain[row.key[1]][0]); \
                if (d < 17) { \
                    e.cons += -100; \
                } else { \
                    e.cons += _mc_constraints_weight * 0.1 * row.weight * square(d - 17); \
                } \
            } \
        } \
    }
    MCPSB_ENERGY_CONSTRAINTS(partial,if (is_selected(row.key[0]) ^ is_selected(row.key[1])))
    MCPSB_ENERGY_CONSTRAINTS(total,)

    double total_energy() {
        en_t e;
        mc_total_energy(e);
        e.print();
    }

};

using EnPsb = En<CGpsb>;

} // namespace mc
} // namespace nuc3d
} // namespace jian

