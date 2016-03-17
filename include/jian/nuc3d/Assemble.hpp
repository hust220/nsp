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

namespace jian {

class Assemble : public BasicPredict3D, public JobPredict3D, public Rand {
public:    
    using Mat = MatrixXd;

    std::map<loop *, std::pair<std::deque<TemplRec>, std::deque<TemplRec>>> _records;
    std::map<loop *, std::pair<Model, Model>> _templates;
    int _it_num = 0;
    std::map<loop *, bool> _is_virtual;
    std::map<loop *, bool> _is_sampled;
    std::map<loop *, bool> _find_self;

    BuildStrand build_strand;
    BuildLoop build_loop;
    BuildLoopDG build_loop_dg;

    Assemble(const Par &par) : JobPredict3D(par) {
        Trace::log("# Check input\n");
        if (!NucSS::check_ss(_ss)) throw "The secondary structure includes illegal characters!";

        Trace::log("# Construct 2D Structure Tree\n");
        _ss_tree.make(_seq, _ss, _hinge);

        Trace::log("# Set Virtual Loops\n");
        set_virtual_loops();

        Trace::log("# Searching Templates\n");
        find_templates();

        Trace::log("# Set the Loops to be Sampled\n");
        set_loops_sampled();

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

    Model predict() {
        _it_num++;
        Trace::log("# Select Templates\n");
        if (_it_num == 1) select_templates(); else sample();
        Trace::log("# Print Templates\n");
        LOOP_TRAVERSE(_ss_tree.head, Trace::log(L, " : ", _templates[L].first.name, ' ', _templates[L].second.name, '\n'));
        Trace::log("# Position Templates\n");
        position_templates();
        Trace::log("# Assemble Templates\n");
        auto &&residues = assemble_templates(_ss_tree.head, _seq.size());
        Trace::log("# Build Strands\n");
        build_strands(residues);
        Trace::log("# Transform\n");
        return Transform(residues_to_model(std::move(residues)))(_type, _seq);
    }

    void select_templates() {
        LOOP_TRAVERSE(_ss_tree.head,
            if (L->has_loop()) {
                if (_records[L].first.empty()) {
                    if (_is_virtual[L]) {
                        _templates[L].first = build_loop(L->seq(), NucSS::hinge_ss(L->ss()));
                    } else {
                        build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); _templates[L].first = build_loop_dg();
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

    void sample() {
        LOOP_TRAVERSE(_ss_tree.head,
            if (L->has_loop() && _is_sampled[L]) {
                int index = int(rand() * _num_sampling);
                if (index >= _records[L].first.size()) {
                    if (_is_virtual[L]) {
                        _templates[L].first = build_loop(L->seq(), NucSS::hinge_ss(L->ss()));
                    } else {
                        build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); _templates[L].first = build_loop_dg();
                    }
                } else {
                    _templates[L].first = load_pdb(_records[L].first[index], (_is_virtual[L] ? "hinge" : ""));
                }
            }
        );
    }

    Model load_pdb(const TemplRec &templ_res, const std::string &type = "") {
        std::string lib_path = _lib + "/" + _type + "/templates/";
        std::string pdb_name = templ_res._name + ((type == "") ? "" : ("." + type)) + ".pdb";
        return Model(lib_path + pdb_name);
    }

    Chain assemble_templates(loop *l, int len) {
        Chain chain; chain.resize(len);
        LOOP_TRAVERSE(l,
            if (L->has_loop()) {
                auto temp_residues = residues(_templates[L].first);
                int i = 0;
                LOOP_EACH(L,
                    if (((!(_is_virtual[L]) && RES->type != '&') || RES->type == '(' || RES->type == '(')) {
                        chain[RES->num - 1] = temp_residues[i];
                        i++;
                    }
                );
            }
            if (L->has_helix()) {
                int len = L->s.len();
                auto temp_residues = residues(_templates[L].second);
                HELIX_EACH(L->s,
                    chain[BP->res1.num - 1] = temp_residues[N_BP];
                    chain[BP->res2.num - 1] = temp_residues[2*len-1-N_BP];
                );
            }
        );
        return chain;
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
        auto c1 = -sp.c1; auto &rot = sp.rot; auto &c2 = sp.c2;
        EACH_ATOM(model, geom::translate(ATOM, c1); geom::rotate(ATOM, rot); geom::translate(ATOM, c2));
    }

    Mat model_mat(const Model &model, const std::list<int> &list) {
        static std::set<std::string> names{"C5*", "O3*", "C1*"};
        Mat mat(names.size() * list.size(), 3);
        int index = 0;
        EACH_RES(model, IF(HAS(list, N_RES), EACH(atom, RES,
            IF(names.count(atom.name), LOOP(i, 3, mat(index, i) = atom[i]); index++))));
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
        Trace::log("make strand ("); for (auto i : vec) Trace::log(i, ' '); Trace::log(")...\n");
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
        Trace::log("Records searching results:\n");
        LOOP_TRAVERSE(_ss_tree.head,
            Trace::log("Hairpin(", L, "):\n", 
                       "  Helix: ", L->s.seq(), ' ', L->s.ss(), ' ', _records[L].second.size(), '\n',
                       "  Loop: ", L->num_sons(), ' ', L->seq(), ' ', L->ss(), ' ', _records[L].first.size(), "\n\n");
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

