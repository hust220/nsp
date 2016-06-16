#pragma once

#include <memory>
#include <sstream>
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../mc.hpp"
#include "../utils/file.hpp"
#include "BuildHelix.hpp"
#include "transform.hpp"
#include "C2A.hpp"
#include "JobPredict3D.hpp"

namespace jian {
namespace nuc3d {

class MC1p : public JobPredict3D, public MC {
public:    
    enum res_module_t {RES_LOOP, RES_HELIX, RES_HAIRPIN};

    struct en_t {
        double len = 0, ang = 0, dih = 0, crash = 0, cons = 0, vdw = 0;
        double sum() const { return len + ang + dih + crash + cons + vdw; }
    };

    using item_t = Atom;
    using space_val_t = std::list<int>;
    using space_t = std::map<int, std::map<int, std::map<int, space_val_t>>>;
    using item_space_t = std::vector<space_val_t *>;
    using range_t = std::array<int, 4>;

    std::deque<std::shared_ptr<range_t>> m_range;
    space_t m_space;
    item_space_t m_item_space;
    int m_box = 2;
    std::vector<res_module_t> _res_module_types;
    std::vector<loop *> _res_module;
    int _mc_index;
    std::deque<std::array<double, 3>> _moved_atoms;
    range_t _moved_residues;
    std::deque<int> _free_atoms;
    std::deque<std::shared_ptr<SSTree>> m_trees;
    std::string m_init_pdb;
    Chain _pred_chain;
    std::vector<std::string> _suppos_atoms {"C4*"};

    double _mc_max_shift = 2;

    double _mc_bond_length_weight = 10;
    double _mc_bond_angle_weight = 5;
    double _mc_bond_dihedral_weight = 3;
    double _mc_constraints_weight = 1;
    double _mc_crash_weight = 100;
    double _mc_vdw_weight = 0.01;

    double _mc_bond_angle_std = 2.6;
    double _mc_bond_dihedral_std = 0.27;

    MC1p(const Par &par) : JobPredict3D(par) {
        par.set(m_init_pdb, "pdb");
        _pred_chain = residues_from_file(m_init_pdb);
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

        std::cout << "# Set MC energy weights." << std::endl;
        par.set(_mc_bond_length_weight, "mc_bond_length_weight");
        par.set(_mc_bond_angle_weight, "mc_bond_angle_weight");
        par.set(_mc_bond_angle_std, "mc_bond_angle_std");
        par.set(_mc_bond_dihedral_weight, "mc_bond_dihedral_weight");
        par.set(_mc_bond_dihedral_std, "mc_bond_dihedral_std");
        par.set(_mc_constraints_weight, "mc_constraints_weight");
        par.set(_mc_crash_weight, "mc_crash_weight");
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

    void set_free_atoms(const Par::val_t &par) {
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
//        _res_module_types.resize(_seq.size(), RES_LOOP);
//        _res_module.resize(_seq.size(), NULL);
//        for (int i = 0; i < _seq.size(); i++) {
//            m_range[i] = std::make_shared<range_t>({i, i, i, i});
//        }
        std::vector<int> v(_seq.size(), 0);

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
//                std::cout << flag << ' ' << n << std::endl;
                if (flag == 0 && n <= 14) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        };

        auto make_range = [](int a, int b, int c,  int d) {
            std::shared_ptr<range_t> p = std::make_shared<range_t>();
            (*p) = {a, b, c, d};
            return p;
        };

        auto make_hairpin_range = [&](auto &&l) {
            return make_range(l->s.head->res1.num - 1, l->s.head->res2.num - 1, l->s.head->res1.num - 1, l->s.head->res2.num - 1 );
        };

        auto make_helix_range = [&](auto &&h) {
            auto p = make_range(h.head->res1.num - 1, 0, 0, h.head->res2.num - 1);
            HELIX_EACH(h,
                if (BP->next == NULL) {
                    (*p)[1] = BP->res1.num - 1;
                    (*p)[2] = BP->res2.num - 1;
                }
            );
            return p;
        };

        auto set_res_module_types_ss = [&](loop *l, bool is_first){
            LOOP_TRAVERSE(l,
                if (!m_sample_hairpin && is_first && is_hairpin(L)) {
                    m_range.push_back(make_hairpin_range(L));
                    for (int i = L->s.head->res1.num - 1; i <= L->s.head->res2.num - 1; i++) {
//                        _res_module_types[i] = RES_HAIRPIN;
//                        _res_module[i] = L;
                        v[i] = 1;
                    }
                } else if (L->has_helix()) {
                    m_range.push_back(make_helix_range(L->s));
                    HELIX_EACH(L->s,
//                        _res_module_types[BP->res1.num - 1] = RES_HELIX;
//                        _res_module[BP->res1.num - 1] = L;
                        v[BP->res1.num - 1] = 1;
//                        _res_module_types[BP->res2.num - 1] = RES_HELIX;
//                        _res_module[BP->res2.num - 1] = L;
                        v[BP->res2.num - 1] = 1;
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

        for (int i = 0; i < _seq.size(); i++) {
            if (v[i] == 0) {
                m_range.push_back(make_range(i, i, i, i));
            }
        }
    }

    void mc_write() {
        static int n = 1;
        std::ostringstream stream;
        stream << _name << ".mc." << m_seed << ".pdb";
        std::string name = stream.str();
        if (n == 1) {
            file::clear_file(name);
        }
        append_chain_to_file(_pred_chain, name, n);
        en_t e;
        mc_total_energy(e);
        std::cout << _mc_step + 1 << ": " <<  e.sum() << "(total) " << e.crash << "(crash) " << e.vdw << "(vdw) " << e.len << "(bond) " << e.ang << "(ang) " << e.dih << "(dih) " << e.cons << "(c) "<< _mc_tempr << "(tempr) " << _mc_local_succ_rate << "(rate)" << std::endl;
        n++;
    }

    void mc_sample() {
        backup();
        int beg = _moved_residues[0];
        int end = _moved_residues[3];
        if (beg == end || rand() < 0.5) {
            // translate
            int index = int(rand() * 3);
            double dist = (rand() - 0.5) * 2 * _mc_max_shift;
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    _pred_chain[i][0][index] += dist;
                    space_update_item(i);
                }
            }
        } else {
            // rotate
            int index = int(rand() * 3);
            double dih = (rand() - 0.5) * PI / 6;
            auto &&rot = geom::rot_mat(index, dih);
            auto &&origin = center(_pred_chain[beg][0], _pred_chain[end][0]);
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    geom::rotate(_pred_chain[i][0], origin, rot);
                    space_update_item(i);
                }
            }
        }
    }

    void mc_back() {
        for (int i = 0; i < _seq.size(); i++) {
            if (is_selected(i)) {
                item_t &atom = item(i);
//                for (auto && atom : _pred_chain[i]) {
                    auto &arr = _moved_atoms.front();
                    atom[0] = arr[0]; atom[1] = arr[1]; atom[2] = arr[2];
                    _moved_atoms.pop_front();
//                }
                space_update_item(i);
            }
        }
    }

    void backup() {
        _moved_atoms.clear();
        for (int i = 0; i < _seq.size(); i++) {
            if (is_selected(i)) {
                for (auto && atom : _pred_chain[i]) {
                    _moved_atoms.push_back({atom[0], atom[1], atom[2]});
                }
            }
        }
    }

    int space_index(double n) const {
        return (n+1000)/5;
    }

    item_t &item(int i) {
        return _pred_chain[i][0];
    }

    space_val_t &space_val(int i) {
        item_t &a = item(i);
        return m_space[space_index(a[0])][space_index(a[1])][space_index(a[2])];
    }

    void space_update_item(int i) {
        space_val_t &n = space_val(i);
//        item_t &a = item(i);
        space_val_t &o = *(m_item_space[i]);
        if (&o != &n) {
            o.erase(std::find(o.begin(), o.end(), i));
            m_item_space[i] = &n;
            n.push_back(i);
        }
    }

    void run() {
        Debug::println("# Set pseudo-knots");
        set_pseudo_knots();
        Debug::print("# Coarse Grained\n");
        _pred_chain = _pred_chain.coarse_grained(_suppos_atoms);
        std::cout << "# Init space..." << std::endl;
        m_item_space.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            space_val_t &s = space_val(i);
//            item_t &a = item(i);
            s.push_back(i);
            m_item_space[i] = &s;
        }
        std::cout << "# MC Start..." << std::endl;
        mc();
        Debug::println("# Print Constraints...");
        print_constraints();
        Debug::println("# Coarsed Grained To All Atom...");
        coarse_grained_to_all_atom();
        std::cout << "# Transform." << std::endl;
        this->transform();
        std::cout << "# Writing to file." << std::endl;
        std::ostringstream stream;
        stream << _name << ".sample." << m_seed << ".pdb";
        residues_to_file(_pred_chain, stream.str());
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
        _pred_chain = c2a(c, 0, num_atoms-1);
    }

    void print_constraints() {
        double d;
        int i, j;
        for (auto && ct : _constraints.distances) {
            i = ct.key[0];
            j = ct.key[1];
            d = geom::distance(_pred_chain[i][0], _pred_chain[j][0]);
            std::cout << i << ' ' << j << " value:" << ct.value << " weight:" << ct.weight << " dist:" << d << std::endl;
        }
    }

    virtual void mc_select() {
        int len = m_range.size();
        _moved_residues = *(m_range[int(rand() * len)]);
    }

    bool is_selected(const int &i) const {
        auto &p = _moved_residues;
        return (i >= p[0] && i <= p[1]) || (i >= p[2] && i <= p[3]);
    }

    double mc_partial_energy() {
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

//        double d; 
//        for (int i = 0; i < _seq.size(); i++) { 
//            for (int j = i + 2; j < _seq.size(); j++) { 
//                cond { 
//                    d = geom::distance(_pred_chain[i][0], _pred_chain[j][0]); 
//                    if (d < 5) { 
//                        e.crash += _mc_crash_weight * square(d - 5); 
//                    } 
//                } 
//            } 
//        } 

#define MC_ENERGY_CRASH(name, cond) \
    void mc_##name##_energy_crash(en_t &e) { \
        double d; \
        int a, b, c; \
        for (int n = 0; n < _seq.size(); n++) { \
            cond { \
                for (int i = -m_box; i <= m_box; i++) { \
                    for (int j = -m_box; j <= m_box; j++) { \
                        for (int k = -m_box; k <= m_box; k++) { \
                            item_t &it = item(n); \
                            a = space_index(it[0])+i; \
                            b = space_index(it[1])+j; \
                            c = space_index(it[2])+k; \
                            space_val_t &s = m_space[a][b][c]; \
                            for (auto && t : s) { \
                                if (n < t) { \
                                    d = geom::distance(item(t), it); \
                                    if (d < 5) { \
                                        e.crash += _mc_crash_weight * square(d - 5); \
                                    } else { \
                                        e.vdw += -1 + _mc_vdw_weight * square(d - 10); \
                                    } \
                                } \
                            } \
                        } \
                    } \
                } \
            } \
        } \
    } 
//MC_ENERGY_CRASH(partial, if (is_selected(i) ^ is_selected(j)))
MC_ENERGY_CRASH(partial, if (is_selected(n)))
MC_ENERGY_CRASH(total, )


#define MC_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d, d1, d2, d3; \
        for (int i = 0; i < _seq.size() - 1; i++) { \
            cond { \
                d = geom::distance(_pred_chain[i][0], _pred_chain[i+1][0]); \
                e.len += _mc_bond_length_weight * square(d - 6.1); \
            } \
        } \
    } 
MC_ENERGY_BOND(partial, if (is_selected(i) ^ is_selected(i+1)));
MC_ENERGY_BOND(total,);

#define MC_ENERGY_ANGLE(name, cond) \
    void mc_##name##_energy_angle(en_t &e) { \
        double d; \
        int len = _seq.size(); \
        for (int i = 0; i < len - 2; i++) { \
            cond { \
                d = geom::angle(_pred_chain[i][0],  \
                                _pred_chain[i+1][0],  \
                                _pred_chain[i+2][0]); \
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
            cond {\
                d = geom::dihedral(_pred_chain[i][0],  \
                                   _pred_chain[i+1][0],  \
                                   _pred_chain[i+2][0],  \
                                   _pred_chain[i+3][0]); \
                d = d - _mc_bond_dihedral_std; \
                d = 3.3 - 4 * std::cos(d) + std::cos(2 * d) - 0.44 * std::cos(3 * d); \
                e.dih += _mc_bond_dihedral_weight * d; \
            }\
        } \
    }
MC_ENERGY_DIHEDRAL(partial, if ((is_selected(i) + is_selected(i+1) + is_selected(i+2) + is_selected(i+3)) % 4 != 0))
MC_ENERGY_DIHEDRAL(total, )


#define MC_ENERGY_CONSTRAINTS(name, cond) \
    void mc_##name##_energy_constraints(en_t &e) { \
        double d; \
        int n = _constraints.distances.size(); \
        for (auto && row : _constraints.distances) { \
            cond { \
                d = geom::distance(_pred_chain[row.key[0]][0], _pred_chain[row.key[1]][0]); \
                if (d > row.value + 3 || d < row.value - 3) { \
                    e.cons += _mc_constraints_weight * row.weight * 9; \
                } else { \
                    e.cons += _mc_constraints_weight * row.weight * square(d - row.value); \
                } \
            } \
        } \
    }
MC_ENERGY_CONSTRAINTS(partial, if (is_selected(row.key[0]) ^ is_selected(row.key[1])))
MC_ENERGY_CONSTRAINTS(total,)

    template<typename T, typename U>
    std::array<double, 3> center(T &&r1, U &&r2) {
        return {(r1[0]+r2[0])/2.0, (r1[1]+r2[1])/2.0, (r1[2]+r2[2])/2.0};
    }

};

} // namespace nuc3d
} // namespace jian

