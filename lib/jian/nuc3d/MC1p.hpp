#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../pp.hpp"
#include "../mc.hpp"
#include "../utils/file.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../cg.hpp"
#include "BuildHelix.hpp"
#include "transform.hpp"
#include "TemplRec.hpp"
#include "JobPredict3D.hpp"

namespace jian {
namespace nuc3d {

class MC1p : public JobPredict3D, public MC {
public:    
    #define MC1p_en_t_m len, ang, dih, crash, cons, vdw
    struct en_t {
        #define MC1p_en_t_m_def(a) double a = 0;
        JN_MAP(MC1p_en_t_m_def, MC1p_en_t_m)
        #define MC1p_en_t_m_sum(a) + a
        double sum() const { return 0 JN_MAP(MC1p_en_t_m_sum, MC1p_en_t_m); }
    };

    using item_t = Atom;
    using space_val_t = std::list<int>;
    using space_t = std::map<int, std::map<int, std::map<int, space_val_t>>>;
    using item_space_t = std::vector<space_val_t *>;
    using range_t = std::array<int, 4>;

    std::vector<int> m_indices;
    space_t m_space;
    item_space_t m_item_space;
    int m_box = 3;
    double m_box_size = 6;
    std::deque<range_t *> m_range;
    std::deque<Atom> _moved_atoms;
    range_t _moved_residues;
    std::deque<std::shared_ptr<SSTree>> m_trees;
    Chain _pred_chain;

    #define JN_MC_PARS1 heat_steps, cool_steps, cycle_steps, write_steps, heat_rate, dec_rate, init_tempr
    #define JN_MC_PARS2 bond_length_weight, bond_angle_weight, bond_angle_std, bond_dihedral_weight, bond_dihedral_std, \
                        constraints_weight, crash_weight, vdw_weight, max_shift

    #define JN_MCPSB_DEF_PAR(a) double PP_CAT(_mc_, a);
    JN_MAP(JN_MCPSB_DEF_PAR, JN_MC_PARS2)

    void read_par() {
        Par temp_par(Env::lib() + "/RNA/pars/nuc3d/mc1p/mc.par");
        #define JN_MCPSB_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_MCPSB_TEMPPAR_SET, JN_MC_PARS2)
    }

    MC1p(const Par &par) : JobPredict3D(par) {
        LOG << "# Read parameters" << std::endl;
        read_par();

        LOG << "# Set parameters" << std::endl;
        #define JN_MCPSB_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_MCPSB_PAR_SET, JN_MC_PARS1, JN_MC_PARS2)

        LOG << "# Print parameters" << std::endl;
        print_parameters();

        LOG << "# Read initial structure" << std::endl;
        chain_read_model(_pred_chain, par.get("pdb"));

        LOG << "# Set indices" << std::endl;
        thread_local static std::map<char, int> m {{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
        m_indices.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            m_indices[i] = m[_seq[i]];
        }

        LOG << "# Set 2D trees" << std::endl;
        set_trees();

        LOG << "# Set ranges" << std::endl;
        set_ranges();

        LOG << "# Print ranges" << std::endl;
        print_ranges();

    }

    ~MC1p() {
        for (auto && r : m_range) {
            delete r;
        }
    }

    void print_ranges() {
        for (auto && range : m_range) {
            for (auto && i : *range) {
                LOG << i << ' ';
            }
            LOG << std::endl;
        }
    }

    void print_parameters() {
        #define JN_MCPSB_TEMP(a) LOG << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;
        JN_MAP(JN_MCPSB_TEMP, JN_MC_PARS1, JN_MC_PARS2)
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

    void set_trees() {
        m_trees.push_back(std::make_shared<SSTree>());
        m_trees.back()->make(_seq, _ss, 2);
        auto & keys = NucSS::instance().paired_keys;
        for (auto it = keys.begin()+1; it != keys.end(); it++) {
            auto ss = partial_ss(_ss, *it);
            if (std::any_of(ss.begin(), ss.end(), [](auto &&c){return c != '.';})) {
                m_trees.push_back(std::make_shared<SSTree>());
                m_trees.back()->make(_seq, ss, 1);
            } else {
                break;
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

        auto it = m_trees.begin();
        LOOP_TRAVERSE((*it)->head(), 
            if (L->has_helix() && L->s.len() == 1) {
                set_pseudo_knots_helix(L->s);
            }
        );
        for (it = m_trees.begin() + 1; it != m_trees.end(); it++) {
            LOOP_TRAVERSE((*it)->head(), 
                if (L->has_helix()) {
                    set_pseudo_knots_helix(L->s);
                }
            );
        }
    }

    auto make_range(int a, int b, int c,  int d) {
//        std::shared_ptr<range_t> p = std::make_shared<range_t>();
        range_t *p = new range_t;
        (*p) = {a, b, c, d};
        return p;
    }

    template<typename T>
    auto make_hp_range(T &&l) {
        return make_range(l->s.head->res1.num - 1, l->s.head->res2.num - 1, l->s.head->res1.num - 1, l->s.head->res2.num - 1 );
    }

    template<typename T>
    auto make_il_range(T &&l) {
        return make_range(l->s.head->res1.num - 1, l->son->s.head->res1.num - 2, l->son->s.head->res2.num, l->s.head->res2.num - 1 );
    }

    template<typename T>
    auto make_helix_range(T &&h) {
        auto p = make_range(h.head->res1.num - 1, 0, 0, h.head->res2.num - 1);
        HELIX_EACH(h,
            if (BP->next == NULL) {
                (*p)[1] = BP->res1.num - 1;
                (*p)[2] = BP->res2.num - 1;
            }
        );
        return p;
    }

    void set_ranges() {
        std::vector<int> v(_seq.size(), 0);

        auto is_hp = [&](loop *l) {
            if (l->has_loop() && l->has_helix() && l->num_sons() == 0) {
                int flag = 0, n = 0;
                LOOP_EACH(l,
                    char &c = _ss[RES->num - 1];
                    if (c != '.'  && c != '(' && c != ')') {
                        flag = 1;
                        break;
                    } else {
                        n++;
                    }
                );
                if (flag == 0 && n <= 14) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        };

        auto is_il = [&](loop *l) {
            if (l->has_loop() && l->has_helix() && l->num_sons() == 1) {
                int flag = 0, n = 0;
                LOOP_EACH(l,
                    char &c = _ss[RES->num - 1];
                    if (c != '.'  && c != '(' && c != ')') {
                        flag = 1;
                        break;
                    } else {
                        n++;
                    }
                );
                if (flag == 0 && n <= 20) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        };

        auto update_range = [&](auto &&range) {
            for (int i = range->at(0); i <= range->at(1); i++) { v[i] = 1; }
            for (int i = range->at(2); i <= range->at(3); i++) { v[i] = 1; }
            m_range.push_back(range);
        };

        auto set_res_module_types_ss = [&](loop *l, bool is_first){
            LOOP_TRAVERSE(l,
                if (!m_sample_hairpin && is_first && is_hp(L)) {
                    update_range(make_hp_range(L));
                } else if (is_first && is_il(L)) {
                    update_range(make_il_range(L));
                } else if (L->has_helix()) {
                    update_range(make_helix_range(L->s));
                }
            );
        };

        for (auto it = m_trees.begin(); it != m_trees.end(); it++) {
            set_res_module_types_ss((*it)->head(), it == m_trees.begin());
        }

        for (int i = 0; i < _seq.size(); i++) {
            if (v[i] == 0) {
                m_range.push_back(make_range(i, i, i, i));
            }
        }

        merge_ranges();
    }

    void merge_ranges() {
        LOG << "# Merge ranges..." << std::endl;
        while (true) {
            int flag = 0;
            for (auto && r1 : m_range) {
                for (auto && r2 : m_range) {
                    auto &range1 = *r1;
                    auto &range2 = *r2;
                    if (r1 != r2) {
                        if (range1[1] < range1[2]) {
                            if (range1[1] + 1 == range2[0] && range2[3] + 1 == range1[2]) {
                                if (range2[1] < range2[2]) {
                                    m_range.push_back(make_range(range1[0], range2[1], range2[2], range1[3]));
                                } else {
                                    m_range.push_back(make_range(range1[0], range1[3], range1[0], range1[3]));
                                }
                                delete r1;
                                delete r2;
                                m_range.erase(std::remove_if(m_range.begin(), m_range.end(), [&](auto && t){return t == r1 || t == r2;}), m_range.end());
                                flag++;
                                break;
                            }
                        }
                    }
                }
                if (flag > 0) {
                    break;
                }
            }
            if (flag == 0) {
                break;
            }
        }
    }

    void mc_write() {
        thread_local static int n = 1;
        std::ostringstream stream;
        stream << _name << ".mc1p." << m_seed << ".traj.pdb";
        std::string name = stream.str();

        double old_len; double old_ang; double old_dih; double old_cons; double old_crash;
        if (n == 1) {
            file::clean(name);
        }
        append_chain_to_file(_pred_chain, name, n);
        en_t e;
        mc_total_energy(e);
        #define MC1p_print(a) << e.a << PP_STRING3((a)) << ' '
        LOG << _mc_step + 1 << ": " <<  e.sum() << "(total) "
                  JN_MAP(MC1p_print, MC1p_en_t_m) << _mc_tempr << "(tempr) " << _mc_local_succ_rate << "(rate)" << std::endl;
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
        display_start_information();

        LOG << "# Set pseudo-knots" << std::endl;
        set_pseudo_knots();

        LOG << "# Coarse Grained" << std::endl;
        _pred_chain = CG1p::chain(_pred_chain);

        LOG << "# Init space" << std::endl;
        init_space();

        LOG << "# MC..." << std::endl;
        mc_heat();
        mc_cool();

        LOG << "# Print Constraints..." << std::endl;
        print_constraints();

        LOG << "# Coarsed Grained To All Atom..." << std::endl;
        coarse_grained_to_all_atom();

        LOG << "# Transform." << std::endl;
        this->transform();

        LOG << "# Writing to file." << std::endl;
        std::ostringstream stream;
        stream << _name << ".mc1p." << m_seed << ".pdb";
        write_chain(_pred_chain, stream.str());

        display_end_information();
    }

    void transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

    void coarse_grained_to_all_atom() {
        int num_atoms = _seq.size();
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
        _pred_chain = CG1p::aa(c, 0, num_atoms-1);
    }

    void print_constraints() {
        double d;
        int i, j;
        for (auto && ct : _constraints.distances) {
            i = ct.key[0];
            j = ct.key[1];
            d = geom::distance(_pred_chain[i][0], _pred_chain[j][0]);
            LOG << i << ' ' << j << " value:" << ct.value << " weight:" << ct.weight << " dist:" << d << std::endl;
        }
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

#define MC_ENERGY_CRASH(name, cond1, cond2) \
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
                                    if (d < 8) { \
                                        e.crash += _mc_crash_weight * square(d - 8); \
                                    } \
                                } \
                            } \
                        } \
                    } \
                } \
            } \
        } \
    }

    MC_ENERGY_CRASH(partial, if (is_selected(n)), if (!(is_selected(t)) && (t - n != 1 && n - t != 1)))
    MC_ENERGY_CRASH(total, , if (t - n > 1))

#define MC_ENERGY_BOND(name, cond) \
    void mc_##name##_energy_bond(en_t &e) { \
        double d; \
        Mat arr = Mat::Zero(3, 3); \
        for (int n = 0; n < _seq.size() - 1; n++) { \
             cond { \
                d = geom::distance(_pred_chain[n][0], _pred_chain[n+1][0]); \
                e.len += _mc_bond_length_weight * square(d - 6.1); \
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
                d = geom::distance(_pred_chain[row.key[0]][0], _pred_chain[row.key[1]][0]); \
                if (d < 17) { \
                    e.cons += _mc_constraints_weight * (-100); \
                } else { \
                    e.cons += _mc_constraints_weight * 0.1 * row.weight * square(d - 17); \
                } \
            } \
        } \
    }
    MC_ENERGY_CONSTRAINTS(partial,if (is_selected(row.key[0]) ^ is_selected(row.key[1])))
    MC_ENERGY_CONSTRAINTS(total,)

    template<typename T, typename U>
    std::array<double, 3> center_residues(T &&r1, U &&r2) {
        std::array<double, 3> origin {0, 0, 0};
        int n_atom = 0;
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

