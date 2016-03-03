#ifndef JIAN_NUC3D_ASSEMBLE_H
#define JIAN_NUC3D_ASSEMBLE_H

#include "../etl.h"
#include "../pdb.h"
#include "../geom.h"
#include "../nuc2d/util.h"
#include "../nuc2d/N2D.h"
#include "BasicPredict3D.h"
#include "TemplRec.h"
#include "Transform.h"
#include "AddPhos.h"
#include "BuildStrand.h"
#include "BuildLoop.h"
#include "JobPredict3D.h"

namespace jian {
namespace nuc3d {

class Assemble : public BasicPredict3D, public JobPredict3D, public Rand {
public:    
    using Mat = MatrixXd;
    using Row = std::tuple<nuc2d::loop *, int, int, int>;

    std::map<nuc2d::loop *, std::pair<std::deque<TemplRec>, std::deque<TemplRec>>> _records;
    std::map<nuc2d::loop *, std::pair<Model, Model>> _templates;
    std::list<Row> _loop_nums_table;
    int _it_num = 0;
    std::map<nuc2d::loop *, bool> _is_virtual;

    BuildStrand build_strand;
    BuildLoop build_loop;

    Assemble(const Par &par) : JobPredict3D(par) {
        build_strand.log = log;

        log("Initialize start...\n",
            "Step 1: Check input.\n");
        if (!nuc2d::check_ss(_ss)) throw "The secondary structure includes illegal characters!";

        log("Step 2: Construct 2D structure tree.\n");
        _n2d.hinge_base_pair_num = _hinge;
        _n2d(_seq, _ss);

        log("Step 3: Searching templates.\n");
        find_templates();

        log("Initialize end...\n");
    }

    Model predict() {
        set_virtual_loops();
        select_templates();
        for (auto &&pair : _templates) log(pair.first, " : ", pair.second.first.name, ' ', pair.second.second.name, '\n');
        position_templates();
        auto residues = assemble_templates(_n2d.head, _seq.size());
        build_strands(residues);
        return Transform(residues_to_model(std::move(residues)))(_type, _seq);
    }

    void set_virtual_loops() {
        _n2d.head->apply([&](nuc2d::loop *l){
            if (l == NULL || l->empty()) return;
            else if (_strategy == "loose" and (l->is_open() or l->num_sons() >= 2)) _is_virtual[l] = true;
            else _is_virtual[l] = false;
        });
    }

    void select_templates() {
        if (_it_num == 0) {
            _n2d.head->apply([&](nuc2d::loop *l){
                if (l->has_loop()) {
                    if (_records[l].first.empty()) {
                        _templates[l].first = build_loop(l->seq(), (_is_virtual[l] ? nuc2d::hinge_ss(l->ss()) : nuc2d::lower_ss(l->ss())));
                    } else {
                        _templates[l].first = load_pdb(_records[l].first[0], (_is_virtual[l] ? "hinge" : ""));    
                    }
                }
                if (l->has_helix()) _templates[l].second = load_pdb(_records[l].second[0]);
            });
        } else {
            _n2d.head->apply([&](nuc2d::loop *l){
                if (l->has_loop()) {
                    if (_records[l].first.empty()) {
                        _templates[l].first = build_loop(l->seq(), (_is_virtual[l] ? nuc2d::hinge_ss(l->ss()) : nuc2d::lower_ss(l->ss())));
                    } else {
                        _templates[l].first = load_pdb(_records[l].first[int(rand()*_records[l].first.size())], (_is_virtual[l] ? "hinge" : ""));
                    }
                }
                if (l->has_helix()) _templates[l].second = load_pdb(_records[l].second[int(rand()*_records[l].second.size())]);
            });
        }
        _it_num++;
    }

    Model load_pdb(const TemplRec &templ_res, const std::string &type = "") {
        std::string lib_path = _lib + "/" + _type + "/templates/";
        std::string pdb_name = templ_res._name + ((type == "") ? "" : ("." + type)) + ".pdb";
        return Model(lib_path + pdb_name);
    }

    std::vector<Residue> assemble_templates(nuc2d::loop *l, int len) {
        std::vector<Residue> residues(len);
        assemble_templates_impl(residues, l);
        return residues;
    }

    void assemble_templates_impl(std::vector<Residue> &residues, nuc2d::loop *l) {
        if (l == NULL) return;
        if (l->has_helix()) {
            int index = 0;
            int len = l->s.size();
            auto temp_residues = jian::residues(_templates[l].second);
            for (nuc2d::bp *p = l->s.head; p != NULL; p = p->next) {
                residues[p->res1.num - 1] = temp_residues[index];
                residues[p->res2.num - 1] = temp_residues[2*len-1-index];
                index++;
            }
        }
        if (l->has_loop()) {
            int index = 0;
            auto temp_residues = jian::residues(_templates[l].first);
            for (nuc2d::res *p = l->head; p != NULL; p = p->next) {
                if (_is_virtual[l] && p->type != '(' && p->type != ')' || p->type == '&') continue;
                residues[p->num - 1] = temp_residues[index];
                index++;
            }
        }
        assemble_templates_impl(residues, l->son);
        assemble_templates_impl(residues, l->brother);
    }

    void position_templates() {
        position_templates(_n2d.head, std::list<Mat>());
    }

    void position_templates(nuc2d::loop *l, std::list<Mat> mats) {
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
        auto c1 = -sp.c1;
        auto &rot = sp.rot;
        auto &c2 = sp.c2;
        for (auto &&chain: model) for (auto &&residue: chain) for (auto &&atom: residue) {
            geom::translate(atom, c1); geom::rotate(atom, rot); geom::translate(atom, c2);
        }
    }

    Mat model_mat(const Model &model, const std::list<int> &list) {
        static std::set<std::string> names{"C5*", "O3*", "C1*"};
        Mat mat(names.size() * list.size(), 3);
        int index = 0;
        int num_res = 0;
        for (auto &&chain: model) {
            for (auto &&residue: chain) {
                if (std::count(list.begin(), list.end(), num_res)) {
                    for (auto &&atom: residue) {
                        if (names.count(atom.name)) {
                            for (int i = 0; i < 3; i++) mat(index, i) = atom[i];
                            index++;
                        }
                    }
                }
                num_res++;
            }
        }
        return mat;
    }

    std::list<Mat> loop_mats(const Model &model, nuc2d::loop *l) {
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

    std::list<std::vector<int>> get_strands(nuc2d::loop *l) {
        std::vector<std::vector<int>> fsm{{1, 0}, {1, 0}};
        auto type_id = [](char c){
            if (std::count_if(std::next(nuc2d::paired_keys.begin(), 1), nuc2d::paired_keys.end(), [&](const std::pair<char, char> &pair){
                return c == pair.first || c == pair.second;
            }) || std::count(nuc2d::unpaired_keys.begin(), nuc2d::unpaired_keys.end(), c)) return 0; 
            else return 1;
        };
        int index = -1, left_index = -1, right_index = -1, state = 0;
        std::list<std::vector<int>> strands;
        l->each([&](nuc2d::res *r, int i){
            if (r->type == '&') return;
            auto new_state = fsm[state][type_id(r->type)];
            index = r->num - 1;
            if (new_state == 1 && state == 0) left_index = index;
            else if (new_state == 0 && state == 1) {
                right_index = index;    
                std::vector<int> vec(right_index - left_index);
                std::iota(vec.begin(), vec.end(), left_index);
                strands.push_back(vec);
            }
            state = new_state;
        });
        if (left_index > right_index) {
            std::vector<int> vec(index - left_index + 1);
            std::iota(vec.begin(), vec.end(), left_index);
            strands.push_back(vec);
        }
        return strands;
    }

    std::vector<Residue> make_strand(const std::vector<Residue> &residues, const std::vector<int> &vec) {
        log("make strand (");
        int len = residues.size();
        for (auto i : vec) log(i, ' ');
        log(")...\n");
        Mat a, b;
        if (vec[0] >= 2) {
            a.resize(2, 3);
            for (int i = 0; i < 3; i++) {
                a(0, i) = atom(residues[vec[0] - 2], "C4*")[i];
                a(1, i) = atom(residues[vec[0] - 1], "C4*")[i];
            }
        }
        if (vec.back() <= len - 3) {
            b.resize(2, 3);
            for (int i = 0; i < 3; i++) {
                b(0, i) = atom(residues[vec.back() + 1], "C4*")[i];
                b(1, i) = atom(residues[vec.back() + 2], "C4*")[i];
            }
        }
        return build_strand(vec.size(), a, b);
    }

    void build_strands(std::vector<Residue> &residues) {
        int len = residues.size();
        _n2d.head->apply([&](nuc2d::loop *l){
            if (_is_virtual[l]) for (auto &&strand : get_strands(l)) {
                auto strand_residues = make_strand(residues, strand);
                for (int i = 0; i < strand.size(); i++) residues[strand[i]] = strand_residues[i];
            }
        });
    }

    void find_templates() {
        find_records(_n2d.head);
        print_records();
    }

    void print_records() {
        log("Records searching results:\n");
        for (auto &&Row: _loop_nums_table) {
            nuc2d::loop *l;
            int type, num_loops, num_helices;
            std::tie(l, type, num_loops, num_helices) = Row;
            log("loop(", l, "):\n", "Helix: ", l->s.seq(), ' ', l->s.ss(), ' ', num_helices, '\n');
            log("Loop: ", l->seq(), ' ', l->ss(), ' ', num_loops, "\n\n");
        }
    }

    void find_records(nuc2d::loop *l) {
        if (l == NULL) return;
        find_loop_records(l);
        find_helix_records(l);
        find_records(l->son);
        find_records(l->brother);
    }

    void find_loop_records(nuc2d::loop *l) {
        if (l->empty()) {
            _loop_nums_table.push_back(std::make_tuple(l, 0, 0, 0));    
            return;
        }

        std::string seq = l->seq(), ss = l->ss(), p_ss = nuc2d::pure_ss(ss), lower_ss = nuc2d::lower_ss(p_ss, 1), family = _family;
        int num_sons = l->num_sons();

        std::string info_file = _lib + "/" + _type + "/" + "records/" + (l->is_open() ? "open_" : "") + "loop";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::nuc3d::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        int num = 0;
        while (ifile >> templ_rec._name >> templ_rec._type >> templ_rec._seq >> templ_rec._ss >> templ_rec._family) {
            if (_strategy == "loose" and templ_rec._type == num_sons and num_sons >= (l->is_open() ? 0 : 2)) {
                templ_rec._score = 0;
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if (not _source_pdb.empty() and templ_rec._name.substr(0, 4) == _source_pdb.substr(0, 4)) {
                    _records[l].first.clear(); _records[l].first.push_back(templ_rec);
                    num = 1;
                    break;
                }
            } else if (nuc2d::pure_ss(nuc2d::lower_ss(templ_rec._ss, 1)) == lower_ss) {
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
                if (not _source_pdb.empty() and templ_rec._name.substr(0, 4) == _source_pdb.substr(0, 4)) {
                    _records[l].first.clear(); _records[l].first.push_back(templ_rec);
                    num = 1;
                    break;
                }
            }
        }
        ifile.close();
        std::sort(_records[l].first.begin(), _records[l].first.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
        _loop_nums_table.push_back(std::make_tuple(l, num_sons, num, 0));
    }

    void find_helix_records(nuc2d::loop *l) {
        if (l->s.empty()) return;

        std::string seq = l->s.seq(), ss = l->s.ss(), family = _family;
        int len = l->s.len();

        std::string info_file = _lib + "/" + _type + "/" + "records/helix";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::nuc3d::FindTemplates error! Can't find '" + info_file + "'!";

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
        for (auto &&tuple: _loop_nums_table) if (std::get<0>(tuple) == l) {std::get<3>(tuple) = num; break;}
    }

};

} // namespace nuc3d
} // namespace jian

#endif

