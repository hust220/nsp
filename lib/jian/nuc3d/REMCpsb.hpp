#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../mc.hpp"
#include "../utils/file.hpp"
#include "../scoring/score_psb.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "BuildHelix.hpp"
#include "transform.hpp"
#include "TemplRec.hpp"
#include "JobPredict3D.hpp"
#include "psb2aa.hpp"
#include "../pp.hpp"

namespace jian {
namespace nuc3d {

class REMCpsb : public JobPredict3D, public REMC {
public:    
//    enum res_module_t {RES_LOOP, RES_HELIX, RES_HAIRPIN};

    struct en_t {
        double len = 0, ang = 0, dih = 0, crash = 0, cons = 0, vdw = 0, stacking = 0, pairing = 0;
        double sum() const { return len + ang + dih + crash + cons + vdw + stacking + pairing; }
    };

    using item_t = Atom;
    using space_val_t = std::list<int>;
    using space_t = std::map<int, std::map<int, std::map<int, space_val_t>>>;
    using item_space_t = std::vector<space_val_t *>;
    using range_t = std::array<int, 4>;

    std::vector<std::vector<double>> m_stacking_pars {
        {5.2, 6.5, 5.3, 4.9, 5.1, 6.7, 5.1, 5.1},
        {4.5, 3.9, 4.3, 3.3, 4.5, 3.9, 4.5, 3.3},
        {5.3, 6.7, 5.1, 4.9, 5.3, 6.7, 5.3, 4.9},
        {5.5, 4.9, 5.3, 3.9, 5.5, 5.1, 5.5, 4.1},
        {5.3, 6.7, 5.3, 5.1, 5.3, 6.9, 5.3, 5.1},
        {4.9, 3.9, 4.9, 3.3, 4.3, 4.1, 4.5, 3.3},
        {5.3, 6.7, 5.3, 5.1, 5.3, 6.7, 5.3, 5.3},
        {5.5, 5.1, 5.7, 4.1, 5.3, 5.1, 5.7, 4.1}
    };

    std::vector<std::vector<double>> m_pairing_pars {
        {15.5, 13.0, 11.4},
        {13.0, 10.8,  9.2},
        { 9.6,  7.4,  5.4},
    };

    std::vector<int> m_indices;
    space_t m_space;
    item_space_t m_item_space;
    int m_box = 3;
    double m_box_size = 6;
    std::deque<std::shared_ptr<range_t>> m_range;
    int _mc_index;
    std::deque<Atom> _moved_atoms;
    range_t _moved_residues;
    std::deque<int> _free_atoms;
    std::deque<std::shared_ptr<SSTree>> m_trees;
    Chain _pred_chain;

    double _mc_max_shift;
    double _mc_bond_length_weight;
    double _mc_bond_angle_weight;
    double _mc_bond_dihedral_weight;
    double _mc_constraints_weight;
    double _mc_crash_weight;
    double _mc_stacking_weight;
    double _mc_pairing_weight;
    double _mc_vdw_weight;
    double _mc_bond_angle_std;
    double _mc_bond_dihedral_std;

    REMCpsb(const Par &par) : JobPredict3D(par) {
        #define JN_REMC_PARS1 heat_steps, cool_steps, cycle_steps, write_steps, heat_rate
        #define JN_REMC_PARS2 bond_length_weight, bond_angle_weight, bond_angle_std, bond_dihedral_weight, bond_dihedral_std, \
                            constraints_weight, crash_weight, pairing_weight, stacking_weight, vdw_weight, max_shift

        Par temp_par(Env::lib() + "/RNA/pars/nuc3d/mcpsb/mc.par");
        #define JN_REMCPSB_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_REMCPSB_TEMPPAR_SET, JN_REMC_PARS2)
        std::cout << PP_STRING3(PP_CAT(mc_, bond_angle_weight)) << std::endl;

        std::map<char, int> m{{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
        _pred_chain = residues_from_file(par["pdb"][0]);

        m_indices.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            m_indices[i] = m[_seq[i]];
        }

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


        #define JN_REMCPSB_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        std::cout << "# Set parameters" << std::endl;
        JN_MAP(JN_REMCPSB_PAR_SET, JN_REMC_PARS1, JN_REMC_PARS2)

        std::cout << "# Print parameters" << std::endl;
        #define JN_REMCPSB_TEMP(a) std::cout << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;
        JN_MAP(JN_REMCPSB_TEMP, JN_REMC_PARS1, JN_REMC_PARS2)

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
                        v[i] = 1;
                    }
                } else if (L->has_helix()) {
                    m_range.push_back(make_helix_range(L->s));
                    HELIX_EACH(L->s,
                        v[BP->res1.num - 1] = 1;
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
        std::cout << _mc_step + 1 << ": " <<  e.sum() << "(total) " << e.crash << "(crash) " << e.len << "(bond) " << e.ang << "(ang) " << e.dih << "(dih) " << e.cons << "(c) "<< e.stacking << "(stacking) "  << e.pairing << "(pairing) " << _mc_tempr << "(tempr) " << _mc_local_succ_rate << "(rate)" << std::endl;
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
                        space_update_item(i);
                    }
                }
            }
        } else {
            // rotate
            int beg = _moved_residues[0];
            int end = _moved_residues[3];
            int index = int(rand() * 3);
            double dih = (rand() - 0.5) * PI / 6;
            auto &&rot = geom::rot_mat(index, dih);
            auto &&origin = center_residues(_pred_chain[beg], _pred_chain[end]);
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    for (auto && atom : _pred_chain[i]) {
                        geom::rotate(atom, origin, rot);
                    }
                    space_update_item(i);
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
                space_update_item(i);
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

    void init_space() {
        m_item_space.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            space_val_t &s = space_val(i);
            s.push_back(i);
            m_item_space[i] = &s;
        }
    }

    int space_index(double n) const {
        return (n+1000)/m_box_size;
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
        _pred_chain = cg_psb_chain(_pred_chain);
        std::cout << "# Init space..." << std::endl;
        init_space();
        Debug::print("# REMC...\n");
        REMC::run();
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
        int num_atoms = _seq.size() * 3;
        Eigen::MatrixXd c(num_atoms, 3);
        int num_atom = 0;
        for (int i = 0; i < _pred_chain.size(); i++) {
            for (auto && atom : _pred_chain[i]) {
                for (int j = 0; j < 3; j++) {
                    c(num_atom, j) = atom[j];
                }
                num_atom++;
            }
        }
        _pred_chain = psb2aa(c, 0, num_atoms-1);
//        std::cout << _pred_chain << std::endl;
    }

    void print_constraints() {
        double d;
        int i, j;
        for (auto && ct : _constraints.distances) {
            i = ct.key[0];
            j = ct.key[1];
            d = geom::distance(residue_center(_pred_chain[i]), residue_center(_pred_chain[j]));
            std::cout << i << ' ' << j << " value:" << ct.value << " weight:" << ct.weight << " dist:" << d << std::endl;
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

    bool is_stacking(const Mat &m, int a, int b) const {
        double d;
        int x = m_indices[a] * 2, y = m_indices[b] * 2;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                d = std::fabs(m(i, j) - m_stacking_pars[x + i][y + j]);
                if (d > 1) {
                    return false;
                }
            }
        }
        return true;
    }

    double en_stacking(const Mat &m, int a, int b) const {
        double en = 0, d;
        int x = m_indices[a] * 2, y = m_indices[b] * 2;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                en += square(m(i, j) - m_stacking_pars[x + i][y + j]);
            }
        }
        return (-10 + en) * _mc_stacking_weight;
    }

    bool is_pairing(const Mat &m, int a, int b) const {
        int x = m_indices[a]+m_indices[b];
        int y = m_indices[a]*m_indices[b];
        double d;
        if (x == 1 || x == 5 || y == 2) {
            if (_seq[a] == 'A' || _seq[a] == 'G') {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        d = std::fabs(m(i, j) - m_pairing_pars[i][j]);
                        if (d > 1) { 
//                            std::cout << a << ' ' << b << ' ' << i << ' ' << j << ' ' << m(i, j) << ' ' << m_pairing_pars[i][j] << ' ' << d << std::endl;
//                            for (int k = 0; k < 3; k++) std::cout << _pred_chain[a][0][k] << ' '; std::cout << std::endl;
//                            for (int k = 0; k < 3; k++) std::cout << _pred_chain[b][0][k] << ' '; std::cout << std::endl;
                            return false;
                        }
                    }
                }
            } else {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        d = std::fabs(m(i, j) - m_pairing_pars[j][i]);
                        if (d > 1.5) {
//                            std::cout << a << ' ' << b << ' ' << i << ' ' << j << ' ' << m(i, j) << ' ' << m_pairing_pars[i][j] << ' ' << d << std::endl;
                            return false;
                        }
                    }
                }
            }
            return true;
        } else {
            return false;
        }
    }

    double en_pairing(const Mat &m, int a, int b) const {
        int x = m_indices[a]+m_indices[b];
        int y = m_indices[a]*m_indices[b];
        double d, en = 0;
        if (x == 1 || x == 5 || y == 2) {
            if (_seq[a] == 'A' || _seq[a] == 'G') {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        en += square(m(i, j) - m_pairing_pars[i][j]);
                    }
                }
            } else {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        en += square(m(i, j) - m_pairing_pars[j][i]);
                    }
                }
            }
            if (x == 1) {
                return (-25+en) * 2 * _mc_pairing_weight;
            } else if (x == 5) {
                return (-25+en) * 3 * _mc_pairing_weight;
            } else {
                return (-25+en) * 0.5 * _mc_pairing_weight;
            }
        } else {
            return 0;
        }
    }

#define MC_ENERGY_CRASH(name, cond) \
    void mc_##name##_energy_crash(en_t &e) { \
        int a, b, c; \
        double d; \
        Mat arr = Mat::Zero(3, 3); \
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
                                if (t - n > 1) { \
                                    for (int i = 0; i < 3; i++) { \
                                        for (int j = 0; j < 3; j++) { \
                                            d = geom::distance(_pred_chain[n][i], _pred_chain[t][j]); \
                                            arr(i, j) = d; \
                                            if ((i == 0 || j == 0) && d < 7) { \
                                                e.crash += _mc_crash_weight * square(d - 7); \
                                            } else if ((i == 1 || j == 1) && d < 5) { \
                                                e.crash += _mc_crash_weight * square(d - 5); \
                                            } else if (d < 3.5) { \
                                                e.crash += _mc_crash_weight * square(d - 3.5); \
                                            } \
                                        } \
                                    } \
                                    if (is_stacking(arr, n, t)) { \
                                        e.stacking += en_stacking(arr, n, t); \
                                    } else if (is_pairing(arr, n, t)) { \
                                        e.pairing += en_pairing(arr, n, t); \
                                    } \
                                } \
                            } \
                        } \
                    } \
                } \
            } \
        } \
    }

    MC_ENERGY_CRASH(partial, if (is_selected(n)))
    MC_ENERGY_CRASH(total, )

#define MC_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d; \
        Mat arr = Mat::Zero(3, 3); \
        for (int n = 0; n < _seq.size() - 1; n++) { \
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
                if (is_stacking(arr, n, n+1)) {  \
                    e.stacking += en_stacking(arr, n, n+1);  \
                } \
             } \
        } \
    }
    MC_ENERGY_BOND(partial, if (is_selected(n) ^ is_selected(n+1)))
    MC_ENERGY_BOND(total, )

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
        int len = m_range.size();
        _moved_residues = *(m_range[int(rand() * len)]);
    }

    virtual bool is_selected(const int &i) const {
        auto &p = _moved_residues;
        return (i >= p[0] && i <= p[1]) || (i >= p[2] && i <= p[3]);
    }

};

} // namespace nuc3d
} // namespace jian

