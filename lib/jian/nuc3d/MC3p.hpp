#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../mc.hpp"
#include "../utils/file.hpp"
#include "../scoring/score_3p.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "BuildHelix.hpp"
#include "transform.hpp"
#include "TemplRec.hpp"
#include "JobPredict3D.hpp"
#include "CG2AA.hpp"


namespace jian {
namespace nuc3d {

class MC3p : public JobPredict3D, public MC {
public:    
    enum res_module_t {RES_LOOP, RES_HELIX, RES_HAIRPIN};
    using en_t = struct {double len = 0, ang = 0, dih = 0, crash = 0, cons = 0, stacking = 0, pairing = 0;};

    std::vector<res_module_t> _res_module_types;
    std::vector<loop *> _res_module;
    int _mc_index;
    std::deque<Atom> _moved_atoms;
    std::vector<int> _moved_residues;
    std::deque<int> _free_atoms;
    std::deque<std::shared_ptr<SSTree>> m_trees;
    Chain _pred_chain;
    std::set<std::string> _suppos_atoms{"C5*", "O3*", "C1*"};

    double _mc_max_shift = 2;

    double _mc_bond_length_weight = 10;
    double _mc_bond_angle_weight = 5;
    double _mc_bond_dihedral_weight = 3;
    double _mc_constraints_weight = 1;
    double _mc_crash_weight = 100;
    double _mc_stacking_weight = 1;
    double _mc_pairing_weight = 5;
    double _mc_vdw_weight = 0.1;

    double _mc_bond_angle_std = 2.6;
    double _mc_bond_dihedral_std = 0.27;

    MC3p(const Par &par) : JobPredict3D(par) {
        chain_read_model(_pred_chain, par.get("pdb"));

        Debug::println("# Set residue module types");
        set_res_module_types();

        Debug::println("# Set free atoms");
        if (par.has("free_atoms")) {
            set_free_atoms(par["free_atoms"]);
        } else {
            for (int i = 0; i < _seq.size(); i++) {
                _free_atoms.push_back(i);
            }
        }

        Debug::println("# Set MC steps");
        par.set(_mc_heat_steps, "mc_heat_steps");
        par.set(_mc_cool_steps, "mc_cool_steps");
        par.set(_mc_cycle_steps, "mc_cycle_steps");
        par.set(_mc_write_steps, "mc_write_steps");

        LOG << "# Set MC energy weights." << std::endl;
        par.set(_mc_bond_length_weight, "mc_bond_length_weight");
        par.set(_mc_bond_angle_weight, "mc_bond_angle_weight");
        par.set(_mc_bond_angle_std, "mc_bond_angle_std");
        par.set(_mc_bond_dihedral_weight, "mc_bond_dihedral_weight");
        par.set(_mc_bond_dihedral_std, "mc_bond_dihedral_std");
        par.set(_mc_constraints_weight, "mc_constraints_weight");
        par.set(_mc_crash_weight, "mc_crash_weight");
        par.set(_mc_pairing_weight, "mc_pairing_weight");
        par.set(_mc_stacking_weight, "mc_stacking_weight");
//        par.set(_mc_highest_rate, "mc_highest_rate");
//        par.set(_mc_lowest_rate, "mc_lowest_rate");
//        par.set(_mc_center_rate, "mc_center_rate");
        par.set(_mc_max_shift, "mc_max_shift");

    }

    template<typename T>
    std::string partial_ss(std::string ss, T &&pair) {
        for (auto && c : ss) {
            if (c == pair.first) {
                c = '(';
            } else if (c == pair.second) {
                c = ')';
            } else {
                c = '.';
            }
        }
        return ss;
    }

    void set_free_atoms(const Par::pars_t &par) {
        std::vector<std::string> v;
        for (auto && s : par) {
            tokenize(s, v, "-");
            if (v.size() == 1) {
                _free_atoms.push_back(std::stoi(v[0]) - 1);
            } else if (v.size() == 2) {
                for (int i = std::stoi(v[0]) - 1; i <= std::stoi(v[1]) - 1; i++) {
                    _free_atoms.push_back(i);
                }
            }
        }
    }

    void set_pseudo_knots() {
        auto set_pseudo_knots_helix = [&](auto &&h){
            auto seq = h.seq();
            Model m = build_helix(seq);
            auto nums = h.nums();
            int i = 0;
            for (auto && n : nums) {
                _pred_chain[n] = std::move(m[0][i]);
                i++;
            }
        };

        auto & keys = NucSS::instance().paired_keys;
        for (auto it = keys.begin()+1; it != keys.end(); it++) {
            auto ss = partial_ss(_ss, *it);
            if (std::any_of(ss.begin(), ss.end(), [](auto &&c){return c != '.';})) {
                SSTree ss_tree; 
                ss_tree.make(_seq, ss);
                LOOP_TRAVERSE(ss_tree.head(), 
                    if (L->has_helix()) {
                        set_pseudo_knots_helix(L->s);
                    }
                );
            } else {
                break;
            }
        }
    }

    void set_res_module_types() {
        _res_module_types.resize(_seq.size(), RES_LOOP);
        _res_module.resize(_seq.size(), NULL);

        auto is_hairpin = [&](loop *l) {
            if (l->has_loop() && l->has_helix()) {
                int flag = 0, n = 0;
                LOOP_EACH(l,
                    char c = _ss[RES->num - 1];
                    if (c != '.'  && c != '(' && c != ')') {
                        flag = 1;
                        break;
                    } else {
                        n++;
                    }
                );
//                LOG << flag << ' ' << n << std::endl;
                if (flag == 0 && n <= 14) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        };

        auto set_res_module_types_ss = [&](loop *l, bool is_first){
            LOOP_TRAVERSE(l,
//                L->print();
//                LOG << "is_hairpin: " << is_hairpin(L) << std::endl;
                if (!m_sample_hairpin && is_first && is_hairpin(L)) {
                    for (int i = L->s.head->res1.num - 1; i <= L->s.head->res2.num - 1; i++) {
                        _res_module_types[i] = RES_HAIRPIN;
                        _res_module[i] = L;
                    }
                } else if (L->has_helix()) {
                    HELIX_EACH(L->s,
                        _res_module_types[BP->res1.num - 1] = RES_HELIX;
                        _res_module[BP->res1.num - 1] = L;
                        _res_module_types[BP->res2.num - 1] = RES_HELIX;
                        _res_module[BP->res2.num - 1] = L;
                    );
                }
            );
        };

        auto & keys = NucSS::instance().paired_keys;
        for (auto it = keys.begin(); it != keys.end(); it++) {
            auto ss = partial_ss(_ss, *it);
            if (std::any_of(ss.begin(), ss.end(), [](auto &&c){return c != '.';})) {
                m_trees.push_back(std::make_shared<SSTree>());
                m_trees.back()->make(_seq, ss);
                set_res_module_types_ss(m_trees.back()->head(), it == keys.begin());
            } else {
                break;
            }
        }
    }

    void mc_write() {
        static int n = 1;
        std::ostringstream stream;
        stream << _name << ".mc." << m_seed << ".pdb";
        std::string name = stream.str();
        if (n == 1) {
            file::clean(name);
        }
        append_chain_to_file(_pred_chain, name, n);
        en_t e;
        mc_total_energy(e);
        LOG << _mc_step + 1 << ": " <<  mc_en_sum(e) << "(e_total) " << e.crash << "(e_crash) " << e.len << "(e_bond) " << e.ang << "(en_ang) " << e.dih << "(en_dih) " << e.cons << "(en_c) "<< e.stacking << "(en_stacking) "  << e.pairing << "(en_pairing) " << _mc_tempr << "(temprature) " << _mc_local_succ_rate << "(success rate)" << std::endl;
        n++;
    }

    void mc_sample() {
        backup();
        if (rand() < 0.5) {
            // translate
            int index = int(rand() * 3);
            double dist = (rand() - 0.5) * 2 * _mc_max_shift;
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    for (auto && atom : _pred_chain[i]) {
                        atom[index] += dist;
                    }
                }
            }
        } else {
            // rotate
            int beg = _moved_residues.front();
            int end = _moved_residues.back();
            int index = int(rand() * 3);
            double dih = (rand() - 0.5) * PI / 6;
            auto &&rot = geom::rot_mat(index, dih);
            auto &&origin = center_residues(_pred_chain[beg], _pred_chain[end]);
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    for (auto && atom : _pred_chain[i]) {
                        geom::rotate(atom, origin, rot);
                    }
                }
            }
        }
    }

    void mc_back() {
        for (int i = 0; i < _seq.size(); i++) {
            if (is_selected(i)) {
                for (auto && atom : _pred_chain[i]) {
                    atom = _moved_atoms.front();
                    _moved_atoms.pop_front();
                }
            }
        }
    }

    void backup() {
        _moved_atoms.clear();
        for (int i = 0; i < _seq.size(); i++) {
            if (is_selected(i)) {
                for (auto && atom : _pred_chain[i]) {
                    _moved_atoms.push_back(atom);
                }
            }
        }
    }

    void run() {
        Debug::println("# Set pseudo-knots");
        set_pseudo_knots();
        Debug::print("# Coarse Grained\n");
        _pred_chain = _pred_chain.coarse_grained(_suppos_atoms);
        Debug::print("# MC...\n");
        mc();
        Debug::println("# Print Constraints...");
        print_constraints();
        Debug::println("# Coarsed Grained To All Atom...");
        coarse_grained_to_all_atom();
        LOG << "# Transform." << std::endl;
        this->transform();
        LOG << "# Writing to file." << std::endl;
        std::ostringstream stream;
        stream << _name << ".sample." << m_seed << ".pdb";
        write_chain(_pred_chain, stream.str());
    }

    void transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

    void coarse_grained_to_all_atom() {
        int num_atoms = _seq.size() * _suppos_atoms.size();
        Eigen::MatrixXd c(num_atoms, 3);
        int num_atom = 0;
        for (int i = 0; i < _pred_chain.size(); i++) {
            _pred_chain[i].sort();
            for (auto && atom : _pred_chain[i]) {
                if (std::find(_suppos_atoms.begin(), _suppos_atoms.end(), atom.name) != _suppos_atoms.end()) {
                    for (int j = 0; j < 3; j++) {
                        c(num_atom, j) = atom[j];
                    }
                    num_atom++;
                }
            }
        }
        _pred_chain = cg2aa(c, 0, num_atoms-1);
    }

    void print_constraints() {
        double d;
        int i, j;
        for (auto && ct : _constraints.distances) {
            i = ct.key[0];
            j = ct.key[1];
            d = geom::distance(residue_center(_pred_chain[i]), residue_center(_pred_chain[j]));
            LOG << i << ' ' << j << " value:" << ct.value << " weight:" << ct.weight << " dist:" << d << std::endl;
        }
    }

    std::array<double, 3> residue_backbone_center(const Residue &r) {
        return {(r[0][0] + r[1][0]) / 2.0, (r[0][1] + r[1][1]) / 2.0, (r[0][2] + r[1][2]) / 2.0};
    }

    template<typename T>
    std::array<double, 3> residue_center(T &&r) {
        std::array<double, 3> arr {0, 0, 0};
        int n_atom = 0;
        for (auto && atom : r) {
            for (int i = 0; i < 3; i++) {
                arr[i] += atom[i];
            }
            n_atom++;
        }
        for (int i = 0; i < 3; i++) arr[i] /= n_atom;
        return arr;
    }

    double mc_en_sum(const en_t &e) {
        return e.len + e.ang + e.dih + e.crash + e.cons + e.stacking + e.pairing;
    }

    double mc_partial_energy() {
        en_t e;
        mc_partial_energy_crash(e);
        mc_partial_energy_bond(e);
        mc_partial_energy_angle(e);
        mc_partial_energy_dihedral(e);
        mc_partial_energy_constraints(e);
        return mc_en_sum(e);
    }

    void mc_total_energy(en_t &e) {
        mc_total_energy_crash(e);
        mc_total_energy_bond(e);
        mc_total_energy_angle(e);
        mc_total_energy_dihedral(e);
        mc_total_energy_constraints(e);
    }

#define MC_ENERGY_CRASH(name, cond) \
    void mc_##name##_energy_crash(en_t &e) { \
        double d1, d2, d3, d4, d5; \
        for (int i = 0; i < _seq.size(); i++) { \
            for (int j = i + 2; j < _seq.size(); j++) { \
                cond { \
                    d1 = geom::distance(_pred_chain[i][0], _pred_chain[j][0]); \
                    if (d1 > 20) continue; \
                    d2 = geom::distance(_pred_chain[i][1], _pred_chain[j][1]); \
                    d3 = geom::distance(_pred_chain[i][2], _pred_chain[j][2]); \
                    d4 = geom::distance(_pred_chain[i][0], _pred_chain[j][1]); \
                    d5 = geom::distance(_pred_chain[i][1], _pred_chain[j][0]); \
                    if (d1 < 9 && d1 > 5 && d2 < 9 && d2 > 5 && d3 < 6.5 && d3 > 5.5) { \
                        /* stacking */ \
                        e.stacking += _mc_stacking_weight * scoring::score_stacking_3p(_pred_chain[i], _pred_chain[j]); \
                    } else if (d1 > 16 && d1 < 20 && d2 > 14 && d2 < 18 && d3 > 10 && d3 < 13 && d4 > 15.7 && d4 < 17.7 && d5 > 15.7 && d5 < 17.7) { \
                        /* pairing */ \
                        e.pairing += _mc_pairing_weight * scoring::score_pairing_3p(_pred_chain[i], _pred_chain[j]); \
                    } else if (d1 < 5 || d2 < 5 || d3 < 5 || d4 < 5 || d5 < 5) { \
                        if (d1 < 5) e.crash += _mc_crash_weight * square(d1 - 5); \
                        if (d2 < 5) e.crash += _mc_crash_weight * square(d2 - 5); \
                        if (d3 < 5) e.crash += _mc_crash_weight * square(d3 - 5); \
                        if (d4 < 5) e.crash += _mc_crash_weight * square(d4 - 5); \
                        if (d5 < 5) e.crash += _mc_crash_weight * square(d5 - 5); \
                    } else { \
                        if (d3 < 11) { \
                            e.crash += _mc_vdw_weight * square(d3 - 11); \
                        } else { \
                            /* e.crash += _mc_vdw_weight * 0.01 * square(d3 - 11); */ \
                        } \
                    } \
                } \
            } \
        } \
    }

    MC_ENERGY_CRASH(partial, if (is_selected(i) ^ is_selected(j)))
    MC_ENERGY_CRASH(total, )

#define MC_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d, d1, d2, d3; \
        for (int i = 0; i < _seq.size() - 1; i++) { \
             cond { \
                d = geom::distance(_pred_chain[i][1], _pred_chain[i+1][0]); \
                d1 = geom::distance(_pred_chain[i][0], _pred_chain[i+1][0]); \
                d2 = geom::distance(_pred_chain[i][1], _pred_chain[i+1][1]); \
                d3 = geom::distance(_pred_chain[i][2], _pred_chain[i+1][2]); \
                e.len += _mc_bond_length_weight * (square(d - 3.2) + square(d1 - 6.1) + square(d2 - 6.1)); \
                if (d3 < 5) { \
                    e.crash += _mc_crash_weight * square(d3 - 5); \
                } else if (d3 < 6.5) { \
                    /* stacking */ \
                    e.stacking += _mc_stacking_weight * scoring::score_stacking_3p(_pred_chain[i], _pred_chain[i+1]); \
                } \
            } \
        } \
    }
    MC_ENERGY_BOND(partial, if (is_selected(i) ^ is_selected(i+1)))
    MC_ENERGY_BOND(total, )

#define MC_ENERGY_ANGLE(name, cond) \
    void mc_##name##_energy_angle(en_t &e) { \
        double d; \
        int len = _seq.size(); \
        for (int i = 0; i < len - 2; i++) { \
             cond { \
                d = geom::angle(residue_backbone_center(_pred_chain[i]),  \
                                residue_backbone_center(_pred_chain[i+1]),  \
                                residue_backbone_center(_pred_chain[i+2])); \
                e.ang += _mc_bond_angle_weight * square(d - _mc_bond_angle_std); \
            } \
        } \
    }
    MC_ENERGY_ANGLE(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2)) % 3 != 0))
    MC_ENERGY_ANGLE(total, )

#define MC_ENERGY_DIHEDRAL(name, cond) \
    void mc_##name##_energy_dihedral(en_t &e) { \
        double d; \
        int len = _seq.size(); \
        for (int i = 0; i < len - 3; i++) { \
            cond { \
                d = geom::dihedral(residue_backbone_center(_pred_chain[i]),  \
                                   residue_backbone_center(_pred_chain[i+1]),  \
                                   residue_backbone_center(_pred_chain[i+2]),  \
                                   residue_backbone_center(_pred_chain[i+3])); \
                d = d - _mc_bond_dihedral_std; \
                d = 3.3 - 4 * std::cos(d) + std::cos(2 * d) - 0.44 * std::cos(3 * d); \
                e.dih += _mc_bond_dihedral_weight * d; \
            } \
        } \
        /* \
        for (auto && row : _constraints.dihedrals) { \
            if ((is_selected(row.key[0]) + is_selected(row.key[1]) + is_selected(row.key[2]) + is_selected(row.key[3])) % 4 != 0) { \
                d = geom::dihedral(residue_center(_pred_chain[row.key[0]]), residue_center(_pred_chain[row.key[1]]), \
                                   residue_center(_pred_chain[row.key[2]]), residue_center(_pred_chain[row.key[3]])); \
                d = 1 - std::cos(d - row.value); \
                e += row.weight * square(d); \
            } \
        } \
        */ \
    }
    MC_ENERGY_DIHEDRAL(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2) + is_selected(i+3)) % 4 != 0))
    MC_ENERGY_DIHEDRAL(total, )

#define MC_ENERGY_CONSTRAINTS(name, cond) \
    void mc_##name##_energy_constraints(en_t &e) { \
        double d; \
        int n = _constraints.distances.size(); \
        for (auto && row : _constraints.distances) { \
            cond { \
                d = geom::distance(residue_center(_pred_chain[row.key[0]]), residue_center(_pred_chain[row.key[1]])); \
                e.cons += _mc_constraints_weight * row.weight * square(d - row.value) / double(n); \
            } \
        } \
    }
    MC_ENERGY_CONSTRAINTS(partial,if (is_selected(row.key[0]) ^ is_selected(row.key[1])))
    MC_ENERGY_CONSTRAINTS(total,)

    template<typename T, typename U>
    std::array<double, 3> center_residues(T &&r1, U &&r2) {
        std::array<double, 3> origin {0, 0, 0};
        int n_atom;
        for (auto && atom : r1) {
            for (int i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        }
        for (auto && atom : r2) {
            for (int i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        }
        for (int i = 0; i < 3; i++) origin[i] /= n_atom;
        return origin;
    }

    virtual void mc_select() {
//        LOG << "mc select." << std::endl;
        int len = _free_atoms.size();
        _mc_index = _free_atoms[int(rand() * len)];
        _moved_residues.resize(4);
        if (_res_module_types[_mc_index] == RES_HAIRPIN) {
            loop *l = _res_module[_mc_index];
            _moved_residues[0] = l->s.head->res1.num - 1;
            _moved_residues[1] = l->s.head->res2.num - 1;
            _moved_residues[2] = l->s.head->res1.num - 1;
            _moved_residues[3] = l->s.head->res2.num - 1;
        } else if (_res_module_types[_mc_index] == RES_HELIX) {
            loop *l = _res_module[_mc_index];
            _moved_residues[0] = l->s.head->res1.num - 1;
            _moved_residues[3] = l->s.head->res2.num - 1;
            HELIX_EACH(l->s,
                if (BP->next == NULL) {
                    _moved_residues[1] = BP->res1.num - 1;
                    _moved_residues[2] = BP->res2.num - 1;
                }
            );
        } else if (_res_module_types[_mc_index] == RES_LOOP) {
            _moved_residues[0] = _mc_index;
            _moved_residues[1] = _mc_index;
            _moved_residues[2] = _mc_index;
            _moved_residues[3] = _mc_index;
        }
    }

    virtual bool is_selected(const int &i) const {
        auto &p = _moved_residues;
        return (i >= p[0] && i <= p[1]) || (i >= p[2] && i <= p[3]);
    }

};

} // namespace nuc3d
} // namespace jian

