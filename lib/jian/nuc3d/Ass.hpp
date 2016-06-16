#pragma once

#include <sstream>
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "TemplRec.hpp"
#include "transform.hpp"
#include "AddPhos.hpp"
#include "BuildStrand.hpp"
#include "BuildLoop.hpp"
#include "BuildLoopDG.hpp"
#include "JobPredict3D.hpp"
#include "CG2AA.hpp"

namespace jian {
namespace nuc3d {

class Assemble : public JobPredict3D {
public:    
    using Mat = Eigen::MatrixXd;

    std::map<loop *, std::pair<std::deque<TemplRec>, std::deque<TemplRec>>> _records;
    std::map<loop *, std::pair<Chain, Chain>> _templates;
    int _it_num = 0;
    std::map<loop *, bool> _is_virtual;
    std::map<loop *, bool> _is_sampled;
    std::map<loop *, bool> _find_self;
    Chain _pred_chain;
    std::set<std::string> _suppos_atoms{"C5*", "O3*", "C1*"};
    bool m_sample = false;

    BuildStrand build_strand;
    BuildLoop build_loop;
    BuildLoopDG build_loop_dg;

    Assemble(const Par &par) : JobPredict3D(par) {
        m_sample = par.has("sample");

        Debug::println("# Check input");
        if (!NucSS::check_ss(_ss)) throw "The secondary structure includes illegal characters!";

        Debug::println("# Construct 2D Structure Tree");
        _ss_tree.make(_seq, _ss, _hinge);

        Debug::println("# Set Virtual Loops");
        set_virtual_loops();

        Debug::println("# Searching Templates");
        find_templates();

        Debug::println("# Set the Loops to be Sampled");
        set_loops_sampled();

    }

    void set_virtual_loops() {
        LOOP_TRAVERSE(_ss_tree.head(), 
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
        LOOP_TRAVERSE(_ss_tree.head(),
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

    void run() {
        std::ostringstream stream;
        int n;

        predict_one();
        std::cout << "# Writing assembly structure." << std::endl;
        stream << _name << ".assemble.pdb";
        residues_to_file(_pred_chain, stream.str());
        if (m_sample) {
            n = 1;
            std::cout << "# Writing sampling structure " << n << std::endl;
            stream.str("");
            stream << _name << ".ass." << n << ".pdb";
            residues_to_file(_pred_chain, stream.str());
            for (n = 2; n <= _num; n++) {
                predict_one();
                std::cout << "# Writing sampling structure " << n << std::endl;
                stream.str("");
                stream << _name << ".ass." << n << ".pdb";
                residues_to_file(_pred_chain, stream.str());
            }
        }
    }

    void predict_one() {
        Debug::print("# Select Templates\n");
        select_templates(); 
        Debug::print("# Print Templates\n");
        LOOP_TRAVERSE(_ss_tree.head(), Debug::print(L, " : ", _templates[L].first.model_name, ' ', _templates[L].second.model_name, '\n'));
        Debug::print("# Position Templates\n");
        position_templates();
        Debug::print("# Assemble Templates\n");
        assemble_templates(_ss_tree.head());
        std::cout << "# Build Strands." << std::endl;
        build_strands(_pred_chain);
        std::cout << "# Transform." << std::endl;
        this->transform();
    }

    void transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

    loop *select_loop() {
        int i = 0;
        int n = int(_templates.size() * rand());
        loop *l;

        for (auto && pair : _templates) {
            if (n == i) {
                l = pair.first;
                if (l->has_loop()) {
                    return l;
                } else {
                    return select_loop();
                }
            }
            i++;
        }
    }

    void select_templates() {
        static int it = 0;

        if (it == 0) {
            LOOP_TRAVERSE(_ss_tree.head(),
                if (L->has_loop()) {
                    if (_records[L].first.empty()) {
    //                    if (_is_virtual[L]) {
    //                        _templates[L].first = build_loop(L->seq(), NucSS::hinge_ss(L->ss()));
    //                    } else {
    //                        build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); 
    //                        _templates[L].first = build_loop_dg();
    //                    }
                        build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); 
                        _templates[L].first = build_loop_dg();
                    } else {
    //                    _templates[L].first = load_pdb(_records[L].first[0], (_is_virtual[L] ? "hinge" : ""));    
                        _templates[L].first = load_pdb(_records[L].first[0]);    
                    }
                }
                if (L->has_helix()) {
                    _templates[L].second = load_pdb(_records[L].second[0]);
                }
            );
        } else {
            loop *l = select_loop();
            if (_records[l].first.empty()) {
                build_loop_dg.init(l->seq(), NucSS::lower_ss(l->ss())); 
                _templates[l].first = build_loop_dg();
            } else {
                int n = int(rand() * _records[l].first.size());
                _templates[l].first = load_pdb(_records[l].first[n]);    
            }
        }

        it++;
    }

    Chain load_pdb(const TemplRec &templ_res, const std::string &type = "") {
        std::string lib_path = _lib + "/" + _type + "/templates/";
        std::string pdb_name = templ_res._name + ((type == "") ? "" : ("." + type)) + ".pdb";
        return residues_from_file(lib_path + pdb_name);
    }

    void assemble_templates(loop *l) {
        _pred_chain.resize(_seq.size());
        LOOP_TRAVERSE(l,
L->print();
            if (L->has_loop()) {
//std::cout << num_residues(_templates[L].first) << std::endl;
                auto &temp_residues = _templates[L].first;
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
                auto &temp_residues = _templates[L].second;
                HELIX_EACH(L->s,
                    _pred_chain[BP->res1.num - 1] = temp_residues[N_BP];
                    _pred_chain[BP->res2.num - 1] = temp_residues[2*len-1-N_BP];
                );
            }
        );
    }

    void position_templates() {
        position_templates(_ss_tree.head(), std::list<Mat>());
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
                int len = helix.size();
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

    void position_model(Chain &chain, const Mat &mat) {
        int len = chain.size();
        auto mat2 = model_mat(chain, std::list<int>{0, 1, len - 2, len - 1});
        auto sp = geom::suppos(mat2, mat);
        INIT_SUPPOS(sp);
        for (auto && res : chain) {
            for (auto && atom : res) {
                APPLY_SUPPOS(atom, sp);
            }
        }
//        EACH_ATOM(model, 
//            APPLY_SUPPOS(ATOM, sp);
//        );
    }

    Mat model_mat(const Chain &chain, const std::list<int> &list) {
        Mat mat(_suppos_atoms.size() * list.size(), 3);
        int index = 0;
        for (int n = 0; n < chain.size(); n++) {
            if (std::find(list.begin(), list.end(), n) != list.end()) {
                for (auto && atom : chain[n]) {
                    if (_suppos_atoms.find(atom.name) != _suppos_atoms.end()) {
                        for (int i = 0; i < 3; i++) {
                            mat(index, i) = atom[i];
                        }
                        index++;
                    }
                }
            }
        }
        return mat;
    }

    std::list<Mat> loop_mats(const Chain &chain, loop *l) {
        std::list<Mat> mats;
        if (l->num_sons() == 0) return mats;
        else if (_is_virtual[l]) {
            if (l->is_open()) {
                for (int i = 0; i < l->num_sons(); i++) {
                    mats.push_back(model_mat(chain, std::list<int>{4*i,4*i+1,4*i+2,4*i+3}));
                }
            } else {
                for (int i = 0; i < l->num_sons(); i++) {
                    mats.push_back(model_mat(chain, std::list<int>{4*i+2,4*i+3,4*i+4,4*i+5}));
                }
            }
        } else {
            for (auto && hinge : l->hinges) {
                int a = hinge.first, b = hinge.second;
                mats.push_back(model_mat(chain, std::list<int>{a-1,a,b,b+1}));
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
        LOOP_TRAVERSE(_ss_tree.head(),
            if (_is_virtual[L]) for (auto &&strand : get_strands(L)) {
                auto strand_residues = make_strand(chain, strand);
                for (int i = 0; i < strand.size(); i++) chain[strand[i]] = strand_residues[i];
            }
        );
    }

    void find_templates() {
        LOOP_TRAVERSE(_ss_tree.head(), 
            _find_self[L] = false;
            find_loop_records(L); 
            find_helix_records(L)
        );
        print_records();
    }

    void print_records() {
        Debug::println("Records searching results:");
        LOOP_TRAVERSE(_ss_tree.head(),
            Debug::println("Hairpin(", L, "):");
            Debug::println("  ", L->s, " (", _records[L].second.size(), ")");
            Debug::println("  ", *L, " (", _records[L].first.size(), ")");
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
//            if (templ_rec._name.substr(0, 4) == _disused_pdb) {
            if (std::find(m_disused_pdbs.begin(), m_disused_pdbs.end(), templ_rec._name.substr(0, 4)) != m_disused_pdbs.end()) {
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

} // namespace nuc3d
} // namespace jian

