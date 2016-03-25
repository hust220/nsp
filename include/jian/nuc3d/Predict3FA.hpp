#pragma once

#include "../etl.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "BasicPredict3D.hpp"
#include "TemplRec.hpp"
#include "Transform.hpp"
#include "AddPhos.hpp"
#include "BuildStrand.hpp"
#include "BuildLoop.hpp"
#include "BuildLoopDG.hpp"
#include "JobPredict3D.hpp"
#include "CG2AA.hpp"

namespace jian {

class Predict3FA : public virtual MC, public BasicPredict3D, public JobPredict3D, public virtual Rand {
public:    
    using Mat = MatrixXd;
    enum res_module_t {RES_LOOP, RES_HELIX};

    std::map<loop *, std::pair<std::deque<TemplRec>, std::deque<TemplRec>>> _records;
    std::map<loop *, std::pair<Model, Model>> _templates;
    int _it_num = 0;
    std::map<loop *, bool> _is_virtual;
    std::map<loop *, bool> _is_sampled;
    std::map<loop *, bool> _find_self;
    std::vector<res_module_t> _res_module_types;
    std::vector<loop *> _res_module;
    Chain _pred_chain;
    std::set<std::string> _suppos_atoms{"C5*", "O3*", "C1*"};
    int _mc_index;
    std::array<int, 2> _moved_residues;
    std::deque<Atom> _moved_atoms;

    BuildStrand build_strand;
    BuildLoop build_loop;
    BuildLoopDG build_loop_dg;

    Predict3FA(const Par &par) : JobPredict3D(par) {
        Debug::println("# Check input");
        if (!NucSS::check_ss(_ss)) throw "The secondary structure includes illegal characters!";

        Debug::println("# Construct 2D Structure Tree");
        _ss_tree.make(_seq, _ss, _hinge);

        Debug::println("# Set Constraints Pseudo-knots");
        if (!(par.has("no_pseudo_knots"))) set_constraints_pseudo_knots();

        Debug::println("# Set Residue Module types");
        set_res_module_types();

        Debug::println("# Set Virtual Loops");
        set_virtual_loops();

        Debug::println("# Searching Templates");
        find_templates();

        Debug::println("# Set the Loops to be Sampled");
        set_loops_sampled();

        Debug::println("# Set MC steps");
        par.set(_mc_heat_steps, "mc_heat_steps");
        par.set(_mc_cool_steps, "mc_cool_steps");
        par.set(_mc_cycle_steps, "mc_cycle_steps");
        par.set(_mc_write_steps, "mc_write_steps");

    }

    void set_constraints_pseudo_knots() {
        auto & keys = NucSS::instance().paired_keys;
        for (auto it = keys.begin()+1; it != keys.end(); it++) {
            auto & pair = *it;
            auto ss = _ss;
            for (auto && c : ss) {
                if (c == pair.first) {
                    c = '(';
                } else if (c == pair.second) {
                    c = ')';
                } else {
                    c = '.';
                }
            }
            Debug::println(ss);
            if (std::any_of(ss.begin(), ss.end(), [](auto &&c){return c != '.';})) {
                SSTree ss_tree; 
                ss_tree.make(_seq, ss);
                LOOP_TRAVERSE(ss_tree.head, 
                    if (L->has_helix()) {
                        set_constraints_helix(L->s);
                    }
                );
            } else {
                break;
            }
        }
    }

    template<typename T>
    void set_constraints_helix(T &&h) {
        int len = 0; 
        std::deque<int> s1, s2; 
        HELIX_EACH(h, 
            len++; 
            s1.push_back(BP->res1.num-1); 
            s2.push_back(BP->res2.num-1)
        );
        for (int i = 0; i < len; i++) for (int j = 0; j < i + 1; j++) {
            _constraints.add_distance(make_distance(s1[i], s1[j], HelixPar::dist_a(j-i)));
            _constraints.add_distance(make_distance(s2[i], s2[j], HelixPar::dist_b(j-i)));
            _constraints.add_distance(make_distance(s1[i], s2[j], HelixPar::dist_c(j-i)));
            _constraints.add_distance(make_distance(s1[j], s2[i], HelixPar::dist_d(j-i)));
            _constraints.add_dihedral(make_dihedral(s1[i], s1[j], s2[j], s2[i], HelixPar::dih(j-i)));
        }
        for (int i = 0; i < len; i++) {
            _constraints.add_distance(make_distance(s1[i], s2[i], HelixPar::instance().dist_bp));
        }
    }

    void set_res_module_types() {
        _res_module_types.resize(_seq.size());
        _res_module.resize(_seq.size());
        LOOP_TRAVERSE(_ss_tree.head,
            if (L->has_helix()) {
                HELIX_EACH(L->s,
                    _res_module_types[BP->res1.num - 1] = RES_HELIX;
                    _res_module[BP->res1.num - 1] = L;
                    _res_module_types[BP->res2.num - 1] = RES_HELIX;
                    _res_module[BP->res2.num - 1] = L;
                );
            }
            if (L->has_loop()) {
                LOOP_EACH(L,
                    if (RES->type != '(' && RES->type != ')' && RES->type != '&') {
                        _res_module_types[RES->num - 1] = RES_LOOP;
                        _res_module[RES->num - 1] = L;
                    }
                );
            }
        );
    }

    void set_virtual_loops() {
        LOOP_TRAVERSE(_ss_tree.head, 
            if (L->has_loop()) {
                if (_method == "FADG") {
                    _is_virtual[L] = false;
//                } else if (_strategy == "loose" && (L->is_open() || L->num_sons() >= 2)) {
//                    _is_virtual[L] = true;
                } else {
                    _is_virtual[L] = false;
                }
            }
        );
    }

    void set_loops_sampled() {
        using X = struct {loop *l; int num_sons; std::size_t num_records;};
        std::deque<X> dq;
        LOOP_TRAVERSE(_ss_tree.head,
            if (L->has_loop()) {
                int num_records = _records[L].first.size();
                if (_source_pdb.empty() || !(_find_self[L])) {
                    dq.push_back({L, L->num_sons(), _records[L].first.size()});
                }
            }
            _is_sampled[L] = false;
        );
        std::sort(dq.begin(), dq.end(), [](auto &&a, auto &&b){return a.num_sons > b.num_sons;});
        if (dq.front().num_sons > 1) {
            for (auto &&i : dq) {
                if (i.num_sons > 1 || i.num_records == 0) {
                    _is_sampled[i.l] = true;
                }
            }
        } else {
            for (auto &&i : dq) {
                _is_sampled[i.l] = true;
            }
        }
    }

    void mc_write() {
        static int n = 1;
        if (n == 1) {
            file::clear_file(_name + ".mc.pdb");
        }
        append_chain_to_file(_pred_chain, _name + ".mc.pdb", n);
        std::cout << energy_constraints() << std::endl;
        Debug::print(_mc_step + 1, ": ",  _mc_en, "(energy) ", _mc_tempr, "(temprature) ", _mc_local_succ_rate, "(success rate)\n");
        n++;
    }

    Model predict() {
        Debug::print("# Select Templates\n");
        select_templates(); 
        Debug::print("# Print Templates\n");
        LOOP_TRAVERSE(_ss_tree.head, Debug::print(L, " : ", _templates[L].first.name, ' ', _templates[L].second.name, '\n'));
        Debug::print("# Position Templates\n");
        position_templates();
        Debug::print("# Assemble Templates\n");
        assemble_templates(_ss_tree.head);
        Debug::print("# Build Strands\n");
        build_strands(_pred_chain);
        Debug::print("# Coarse Grained\n");
        _pred_chain = _pred_chain.coarse_grained(_suppos_atoms);
        Debug::print("# MC Heat...\n");
        mc_heat();
        Debug::print("# MC Cool...\n");
        mc_cool();
        Debug::println("# Print Constraints...");
        print_constraints();
        Debug::println("# Coarsed Grained To All Atom...");
        coarse_grained_to_all_atom();
        return Model{};
    }

    void coarse_grained_to_all_atom() {
        int num_atoms = _seq.size() * _suppos_atoms.size();
        MatrixXd c(num_atoms, 3);
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
        auto && chain = CG2AA::cg2aa(c, std::array<int, 2>{0, num_atoms-1});
        residues_to_file(chain, _name + ".final.pdb");
    }

    void print_constraints() {
        double d;
        int i, j;
        for (auto && ct : _constraints.distances) {
            i = ct.key[0];
            j = ct.key[1];
            d = geom::distance(residue_center(_pred_chain[i]), residue_center(_pred_chain[j]));
            std::cout << i << ' ' << j << ' ' << ct.value << ' ' << d << std::endl;
        }
    }

    void select_templates() {
        LOOP_TRAVERSE(_ss_tree.head,
            if (L->has_loop()) {
                if (_records[L].first.empty()) {
                    if (_is_virtual[L]) {
                        _templates[L].first = build_loop(L->seq(), NucSS::hinge_ss(L->ss()));
                    } else {
                        build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); 
                        _templates[L].first = build_loop_dg();
                    }
                } else {
                    _templates[L].first = load_pdb(_records[L].first[0], (_is_virtual[L] ? "hinge" : ""));    
                }
            }
            if (L->has_helix()) {
                _templates[L].second = load_pdb(_records[L].second[0]);
            }
        );
    }

    void mc_select() {
        int len = _seq.size();
        _mc_index = int(rand() * len);
        if (_res_module_types[_mc_index] == RES_HELIX) {
            loop *l = _res_module[_mc_index];
            if (l == _ss_tree.head) {
                mc_select();
                return;
            }
            int beg = l->s.head->res1.num - 1;
            int end = l->s.head->res2.num - 1;
            _moved_residues[0] = beg;
            _moved_residues[1] = end;
        } else if (_res_module_types[_mc_index] == RES_LOOP) {
            _moved_residues[0] = _mc_index;
            _moved_residues[1] = _mc_index;
        }
    }

    double mc_energy() {
        return rand() * 20;
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
        int len = _seq.size();
        double e = 0;
        double d;
        int &beg = _moved_residues[0];
        int &end = _moved_residues[1];
        for (int i = beg; i <= end; i++) {
            for (int j = 0; j <= beg - 1; j++) {
                d = geom::distance(residue_center(_pred_chain[i]), residue_center(_pred_chain[j]));
                if (i - j == 1 && d < 4) e += 9999999;
                else if (i - j == 1 && d > 8) e += 9999999;
                else if (i - j == 2) e += 10 * square(d - 11.7);
                else if (d < 9) e += 9999999;
            }
            for (int j = end + 1; j < len; j++) {
                d = geom::distance(residue_center(_pred_chain[i]), residue_center(_pred_chain[j]));
                if (j - i == 1 && d < 4) e += 9999999;
                else if (j - i == 1 && d > 8) e += 9999999;
                else if (j - i == 2) e += 10 * square(d - 11.7);
                else if (d < 9) e += 9999999;
            }
        }
        if (beg > 0) {
            d = geom::distance(_pred_chain[beg-1]["O3*"], _pred_chain[beg]["C5*"]);
            e += 100 * square(d - 3.0);
        }
        if (end < _seq.size()-1) {
            d = geom::distance(_pred_chain[end]["O3*"], _pred_chain[end+1]["C5*"]);
            e += 100 * square(d - 3.0);
        }
        mc_partial_energy_distance(e, d, beg, end);
        mc_partial_energy_dihedral(e, d, beg, end);
        return e;
    }

    void mc_partial_energy_distance(double &e, double &d, int &beg, int &end) {
        for (auto && row : _constraints.distances) {
            if ((row.key[0] >= beg && row.key[0] <= end) ^
                (row.key[1] >= beg && row.key[1] <= end)) {
                d = geom::distance(residue_center(_pred_chain[row.key[0]]), residue_center(_pred_chain[row.key[1]]));
                e += square(d - row.value);
            }
        }
    }

    void mc_partial_energy_dihedral(double &e, double &d, int &beg, int &end) {
        for (auto && row : _constraints.dihedrals) {
            if (((row.key[0] >= beg && row.key[0] <= end) +
                 (row.key[1] >= beg && row.key[1] <= end) +
                 (row.key[2] >= beg && row.key[2] <= end) + 
                 (row.key[3] >= beg && row.key[3] <= end)) % 4 != 0) {
                d = geom::dihedral(residue_center(_pred_chain[row.key[0]]), residue_center(_pred_chain[row.key[1]]),
                                   residue_center(_pred_chain[row.key[2]]), residue_center(_pred_chain[row.key[3]]));
                d = d - row.value;
                shrink_dihedral(d);
                e += 10 * square(d);
            }
        }
    }

    void shrink_dihedral(double &n) {
        while (true) {
            if (n >= PI) n -= 2 * PI;
            else if (n < -PI) n += 2 * PI;
            else if (n < 0) n = -n;
            else break;
        }
    }

    double energy_constraints() {
        double e = 0;
        double d;
        for (auto && row : _constraints.distances) {
            d = geom::distance(residue_center(_pred_chain[row.key[0]]), residue_center(_pred_chain[row.key[1]]));
            e += square(d - row.value);
        }
        return e;
    }

    void mc_sample() {
        backup();
        int beg = _moved_residues[0];
        int end = _moved_residues[1];
        if (rand() < 0.5) {
            // translate
            int index = int(rand() * 3);
            double dist = (rand() - 0.5) * 2;
            for (int i = beg; i <= end; i++) {
                for (auto && atom : _pred_chain[i]) {
                    atom[index] += dist;
                }
            }
        } else {
            // rotate
            int index = int(rand() * 3);
            double dih = (rand() - 0.5) * PI / 6;
            auto &&rot = geom::rot_mat(index, dih);
            auto &&origin = center_residues(_pred_chain[beg], _pred_chain[end]);
            for (int i = beg; i <= end; i++) {
                for (auto && atom : _pred_chain[i]) {
                    geom::rotate(atom, origin, rot);
                }
            }
        }
    }

    void mc_back() {
        for (int i = _moved_residues[0]; i <= _moved_residues[1]; i++) {
            for (auto && atom : _pred_chain[i]) {
                atom = _moved_atoms.front();
                _moved_atoms.pop_front();
            }
        }
    }

    void backup() {
        _moved_atoms.clear();
        for (int i = _moved_residues[0]; i <= _moved_residues[1]; i++) {
            for (auto && atom : _pred_chain[i]) {
                _moved_atoms.push_back(atom);
            }
        }
    }

    template<typename T, typename U>
    std::array<double, 3> center_residues(T &&r1, U &&r2) {
        std::array<double, 3> origin {0, 0, 0};
        int n_atom;
        EACH(atom, r1,
            for (int i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        );
        EACH(atom, r2,
            for (int i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        );
        for (int i = 0; i < 3; i++) origin[i] /= n_atom;
        return origin;
    }

    Model load_pdb(const TemplRec &templ_res, const std::string &type = "") {
        std::string lib_path = _lib + "/" + _type + "/templates/";
        std::string pdb_name = templ_res._name + ((type == "") ? "" : ("." + type)) + ".pdb";
        return Model(lib_path + pdb_name);
    }

    void assemble_templates(loop *l) {
        _pred_chain.resize(_seq.size());
        LOOP_TRAVERSE(l,
            if (L->has_loop()) {
                auto temp_residues = residues(_templates[L].first);
                int i = 0;
                LOOP_EACH(L,
                    if (((!(_is_virtual[L]) && RES->type != '&') || RES->type == '(' || RES->type == '(')) {
                        _pred_chain[RES->num - 1] = temp_residues[i];
                        i++;
                    }
                );
            }
            if (L->has_helix()) {
                int len = L->s.len();
                auto temp_residues = residues(_templates[L].second);
                HELIX_EACH(L->s,
                    _pred_chain[BP->res1.num - 1] = temp_residues[N_BP];
                    _pred_chain[BP->res2.num - 1] = temp_residues[2*len-1-N_BP];
                );
            }
        );
    }

    void position_templates() {
        position_templates(_ss_tree.head, std::list<Mat>());
    }

    void position_templates(loop *l, std::list<Mat> mats) {
        if (l == NULL) return;
        auto &loop = _templates[l].first;
        auto &helix = _templates[l].second;
        if (l->has_helix()) {
            // position helix
            if (!mats.empty()) {
                auto mat1 = mats.front();
                position_model(helix, mat1);
            }
            if (l->has_loop()) {
                // position loop
                int len = num_residues(helix);
                auto mat2 = model_mat(helix, std::list<int>{len/2-2, len/2-1, len/2, len/2+1});
                position_model(loop, mat2);
                // position son
                if (l->has_son()) {
                    auto son_mats = loop_mats(loop, l);
                    position_templates(l->son, son_mats);
                }
                // position brother
                if (l->has_brother()) {
                    mats.pop_front();
                    position_templates(l->brother, mats);
                }
            } else {
                // position brother
                if (l->has_brother()) {
                    mats.pop_front();
                    position_templates(l->brother, mats);
                }
            }
        } else {
            // position son
            if (l->has_son()) {
                auto son_mats = loop_mats(loop, l);
                position_templates(l->son, son_mats);
            }
        }
    }

    void position_model(Model &model, const Mat &mat) {
        int len = num_residues(model);
        auto mat2 = model_mat(model, std::list<int>{0, 1, len - 2, len - 1});
        auto sp = geom::suppos(mat2, mat);
        INIT_SUPPOS(sp);
        EACH_ATOM(model, 
            APPLY_SUPPOS(ATOM, sp);
        );
    }

    Mat model_mat(const Model &model, const std::list<int> &list) {
        Mat mat(_suppos_atoms.size() * list.size(), 3);
        int index = 0;
        EACH_RES(model, IF(HAS(list, N_RES), EACH(atom, RES,
            if (_suppos_atoms.find(atom.name) != _suppos_atoms.end()) {
                for (int i = 0; i < 3; i++) {
                    mat(index, i) = atom[i];
                }
                index++;
            }
        )));
        return mat;
    }

    std::list<Mat> loop_mats(const Model &model, loop *l) {
        std::list<Mat> mats;
        if (l->num_sons() == 0) return mats;
        else if (_is_virtual[l]) {
            if (l->is_open()) {
                for (int i = 0; i < l->num_sons(); i++) {
                    mats.push_back(model_mat(model, std::list<int>{4*i,4*i+1,4*i+2,4*i+3}));
                }
            } else {
                for (int i = 0; i < l->num_sons(); i++) {
                    mats.push_back(model_mat(model, std::list<int>{4*i+2,4*i+3,4*i+4,4*i+5}));
                }
            }
        } else {
            for (auto && hinge : l->hinges) {
                int a = hinge.first, b = hinge.second;
                mats.push_back(model_mat(model, std::list<int>{a-1,a,b,b+1}));
            }
    //            int a = l->hinges[0].first, b = l->hinges[0].second;
    //            mats.push_back(model_mat(model, std::list<int>{a-1,a,b,b+1}));
        }
        return mats;
    }

    std::list<std::vector<int>> get_strands(loop *l) {
        std::vector<std::vector<int>> fsm{{1, 0}, {1, 0}};
        auto type_id = [](char c){
            if (std::count_if(std::next(NucSS::instance().paired_keys.begin(), 1), NucSS::instance().paired_keys.end(), [&](const std::pair<char, char> &pair){
                return c == pair.first || c == pair.second;
            }) || std::count(NucSS::instance().unpaired_keys.begin(), NucSS::instance().unpaired_keys.end(), c)) return 0; 
            else return 1;
        };
        int index = -1, left_index = -1, right_index = -1, state = 0;
        std::list<std::vector<int>> strands;
        LOOP_EACH(l,
            if (RES->type != '&') {
                auto new_state = fsm[state][type_id(RES->type)];
                index = RES->num - 1;
                if (new_state == 1 && state == 0) left_index = index;
                else if (new_state == 0 && state == 1) {
                    right_index = index;    
                    std::vector<int> vec(right_index - left_index);
                    std::iota(vec.begin(), vec.end(), left_index);
                    strands.push_back(vec);
                }
                state = new_state;
            }
        );
        if (left_index > right_index) {
            std::vector<int> vec(index - left_index + 1);
            std::iota(vec.begin(), vec.end(), left_index);
            strands.push_back(vec);
        }
        return strands;
    }

    template<typename T, typename U>
    Chain make_strand(T &&residues, U &&vec) {
        int len = residues.size();
        Debug::print("make strand ("); for (auto i : vec) Debug::print(i, ' '); Debug::println(")...");
        Mat a, b;
        if (vec[0] >= 2) {
            a.resize(2, 3);
            for (int i = 0; i < 3; i++) {
                a(0, i) = residues[vec[0] - 2]["C4*"][i];
                a(1, i) = residues[vec[0] - 1]["C4*"][i];
            }
        }
        if (vec.back() <= len - 3) {
            b.resize(2, 3);
            for (int i = 0; i < 3; i++) {
                b(0, i) = residues[vec.back() + 1]["C4*"][i];
                b(1, i) = residues[vec.back() + 2]["C4*"][i];
            }
        }
        return build_strand(vec.size(), a, b);
    }

    void build_strands(Chain &chain) {
        int len = chain.size();
        LOOP_TRAVERSE(_ss_tree.head,
            if (_is_virtual[L]) for (auto &&strand : get_strands(L)) {
                auto strand_residues = make_strand(chain, strand);
                for (int i = 0; i < strand.size(); i++) chain[strand[i]] = strand_residues[i];
            }
        );
    }

    void find_templates() {
        LOOP_TRAVERSE(_ss_tree.head, 
            _find_self[L] = false;
            find_loop_records(L); 
            find_helix_records(L)
        );
        print_records();
    }

    void print_records() {
        Debug::println("Records searching results:");
        LOOP_TRAVERSE(_ss_tree.head,
            Debug::println("Hairpin(", L, "):\n", 
                       "  Helix: ", L->s.seq(), ' ', L->s.ss(), ' ', _records[L].second.size(), '\n',
                       "  Loop: ", L->num_sons(), ' ', L->seq(), ' ', L->ss(), ' ', _records[L].first.size(), "\n");
        );
    }

    void find_loop_records(loop *l) {
        if (l->empty()) return; else if (_method == "FADG") return;

        std::string seq = l->seq(), ss = l->ss(), p_ss = NucSS::pure_ss(ss), lower_ss = NucSS::lower_ss(p_ss, 1), family = _family;
        int num_sons = l->num_sons();

        std::string info_file = _lib + "/" + _type + "/" + "records/" + (l->is_open() ? "open_" : "") + "loop";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        int num = 0;
        while (ifile >> templ_rec._name >> templ_rec._type >> templ_rec._seq >> templ_rec._ss >> templ_rec._family) {
            if (templ_rec._name.substr(0, 4) == _disused_pdb) {
                continue;
            } else if (_is_virtual[l] && templ_rec._type == num_sons) {
                templ_rec._score = 0;
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if ((!(_source_pdb.empty())) && (templ_rec._name.substr(0, 4) == _source_pdb.substr(0, 4))) {
                    _records[l].first.clear(); 
                    _records[l].first.push_back(templ_rec);
                    num = 1;
                    _find_self[l] = true;
                    break;
                }
            } else if (NucSS::pure_ss(NucSS::lower_ss(templ_rec._ss, 1)) == lower_ss) {
                templ_rec._score = (templ_rec._ss == ss ? 5 : 0);
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    if (_is_test) continue; else templ_rec._score += 2;
                }
                for (int i = 0; i < templ_rec._seq.size(); i++) {
                    if (seq[i] == templ_rec._seq[i]) templ_rec._score++;
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if ((!(_source_pdb.empty())) && (templ_rec._name.substr(0, 4) == _source_pdb.substr(0, 4))) {
                    _records[l].first.clear(); 
                    _records[l].first.push_back(templ_rec);
                    num = 1;
                    _find_self[l] = true;
                    break;
                }
            }
        }
        ifile.close();
        std::sort(_records[l].first.begin(), _records[l].first.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
//        _loop_nums_table.push_back(std::make_tuple(l, num_sons, num, 0));
    }

    void find_helix_records(loop *l) {
        if (l->s.empty()) return;

        std::string seq = l->s.seq(), ss = l->s.ss(), family = _family;
        int len = l->s.len();

        std::string info_file = _lib + "/" + _type + "/" + "records/helix";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        int num = 0;
        while (ifile >> templ_rec._name >> templ_rec._len >> templ_rec._seq >> templ_rec._ss >> templ_rec._family) {
            if (templ_rec._len == len) {
                templ_rec._score = 0;
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    if (_is_test) continue; else templ_rec._score += 2;
                }
                for (int i = 0; i < templ_rec._seq.size(); i++) {
                    if (seq[i] == templ_rec._seq[i]) {
                        templ_rec._score++;
                    }
                }
                _records[l].second.push_back(std::move(templ_rec));
                num++;
            }
        }
        ifile.close();
        std::sort(_records[l].second.begin(), _records[l].second.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
//        for (auto &&tuple: _loop_nums_table) if (std::get<0>(tuple) == l) {std::get<3>(tuple) = num; break;}
    }

};

} // namespace jian

