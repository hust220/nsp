#pragma once

#include <iostream>
#include <set>
#include <sstream>
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../utils/rand.hpp"
#include "TemplRec.hpp"
#include "transform.hpp"
#include "BuildLoopDG.hpp"
#include "JobPredict3D.hpp"

namespace jian {
namespace nuc3d {

class Assemble : public JobPredict3D {
public:    
    using Mat = Eigen::MatrixXd;

    std::map<loop *, std::pair<std::deque<TemplRec>, std::deque<TemplRec>>> _records;
    std::map<loop *, std::pair<Chain, Chain>> _templates;
    std::map<loop *, std::pair<TemplRec, TemplRec>> m_selected_record;
    int _it_num = 0;
    std::map<loop *, bool> _is_virtual;
    std::map<loop *, bool> _is_sampled;
    std::map<loop *, bool> _find_self;
    Chain _pred_chain;
    std::set<std::string> _suppos_atoms{"C5*", "O3*", "C1*"};
    bool m_sample = false;
    std::string m_sample_mode = "sample_one";

    BuildLoopDG build_loop_dg;

    Assemble(const Par &par) : JobPredict3D(par) {
        m_sample = par.has("sample");
        par.set(m_sample_mode, "sample_mode");

        LOG << "# Check input" << std::endl;
        if (!NucSS::check_ss(_ss)) throw "The secondary structure includes illegal characters!";

        LOG << "# Construct 2D Structure Tree" << std::endl;
        _ss_tree.make(_seq, _ss, _hinge);

        LOG << "# Set Virtual Loops" << std::endl;
        set_virtual_loops();

        LOG << "# Searching Templates..." << std::endl;
        find_templates();

        LOG << "# Printing records..." << std::endl;
        print_records();

        LOG << "# Set the Loops to be Sampled" << std::endl;
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
        LOG << "# Writing assembly structure." << std::endl;
        if (_par->has("out")) {
            residues_to_file(_pred_chain, (*_par)["out"][0]);
        } else {
            stream << _name << ".assemble.pdb";
            residues_to_file(_pred_chain, stream.str());
        }
        if (m_sample) {
            n = 1;
            LOG << "# Writing sampling structure " << n << std::endl;
            stream.str("");
            stream << _name << ".assemble." << n << ".pdb";
            residues_to_file(_pred_chain, stream.str());
            for (n = 2; n <= _num; n++) {
                if (m_sample_mode == "sample_one") {
                    sample_one_template();
                } else if (m_sample_mode == "sample_all") {
                    sample_all_templates();
                } else {
                    throw "Illegal sampling mode!";
                }
                assemble();
                LOG << "# Writing sampling structure " << n << std::endl;
                stream.str("");
                stream << _name << ".assemble." << n << ".pdb";
                residues_to_file(_pred_chain, stream.str());
            }
        }
    }

    void predict_one() {
        LOG << "# Select Templates" << std::endl;
        select_templates(); 

        LOG << "# Print Templates" << std::endl;
        print_templates();

        LOG << "# Assemble" << std::endl;
        assemble();
    }

    double score_templates() {
        double e = 0;
        for (auto && pair : _records) {
            auto & l = pair.first;
            if (l->has_loop()) {
                e += m_selected_record[l].first._score;
            }
            if (l->has_helix()) {
                e += m_selected_record[l].second._score;
            }
        }
        return e;
    }

    bool lack_templates() {
        loop *l;
        for (auto && pair : _templates) {
            l = pair.first;
            if (l->has_loop()) {
                if (_records[l].first.empty()) return true;
            }
        }
        return false;
    }

    void print_templates() {
        LOOP_TRAVERSE(_ss_tree.head(), LOG << L << " : " << _templates[L].first.model_name << ' ' << _templates[L].second.model_name << std::endl;);
    }

    void assemble() {
        LOG << "## Score of Templates: " << score_templates() << std::endl;

        LOG << "## Print Templates." << std::endl;
        print_templates();

        LOG << "## Position Templates." << std::endl;
        position_templates();

        LOG << "## Assemble Templates." << std::endl;
        assemble_templates(_ss_tree.head());

        LOG << "## Transform." << std::endl;
        this->transform();
    }

    void transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

    void select_templates() {
        LOOP_TRAVERSE(_ss_tree.head(),
            if (L->has_loop()) {
                if (_records[L].first.empty()) {
                    build_loop_dg.init(L->seq(), NucSS::lower_ss(L->ss())); 
                    _templates[L].first = build_loop_dg();
                    m_selected_record[L].first = TemplRec{};
                } else {
                    _templates[L].first = load_pdb(_records[L].first[0]);    
                    m_selected_record[L].first = _records[L].first[0];
                }
            }
            if (L->has_helix()) {
                _templates[L].second = load_pdb(_records[L].second[0]);
                m_selected_record[L].second = _records[L].second[0];
            }
        );
    }

    void sample_one_template() {
        sample_loop_template(select_loop());
    }

    void sample_all_templates() {
        for (auto && pair : _templates) {
            loop *l = pair.first;
            if (l->has_loop()) {
                sample_loop_template(l);
            }
        }
    }

    void sample_loop_template(loop *l) {
        if (_records[l].first.empty()) {
            build_loop_dg.init(l->seq(), NucSS::lower_ss(l->ss())); 
            _templates[l].first = build_loop_dg();
            m_selected_record[l].first = TemplRec{};
        } else {
            int n = int(rand() * _records[l].first.size());
            _templates[l].first = load_pdb(_records[l].first[n]);    
            m_selected_record[l].first = _records[l].first[n];
        }
    }

    loop *select_loop() {
        std::deque<loop *> ls;
        for (auto && pair : _records) {
            auto &l = pair.first;
            if (l->has_loop()) {
                ls.push_back(l);
            }
        }
        assert(!ls.empty());
        int n = int(ls.size() * jian::rand());
        return ls[n];
    }

    Chain load_pdb(const TemplRec &templ_res, const std::string &type = "") {
        std::string lib_path = _lib + "/RNA/templates/";
        std::string pdb_name = templ_res._name + ((type == "") ? "" : ("." + type)) + ".pdb";
        return residues_from_file(lib_path + pdb_name);
    }

    void assemble_templates(loop *l) {
        _pred_chain.resize(_seq.size());
        LOOP_TRAVERSE(l,
//L->print();
            if (L->has_loop()) {
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

    void find_templates() {
        LOOP_TRAVERSE(_ss_tree.head(), 
            _find_self[L] = false;
            find_loop_records(L); 
            find_helix_records(L)
        );
    }

    void print_records() {
        LOG << "Records searching results:" << std::endl;
        LOOP_TRAVERSE(_ss_tree.head(),
            LOG << "Hairpin(" << L << "):" << std::endl;
            if (L->has_helix()) {
                LOG << "  Helix " << std::string(L->s) << " (" << _records[L].second.size() << ")" << std::endl;
            }
            if (L->has_loop()) {
                LOG << "  Loop " << std::string(*L) << " (" << _records[L].first.size() << ")" << std::endl;
            }
        );
    }

    void find_loop_records(loop *l) {
        if (l->empty()) return; else if (_method == "FADG") return;

        std::string seq = l->seq(), ss = l->ss(), p_ss = NucSS::pure_ss(ss), lower_ss = NucSS::lower_ss(p_ss, 1), family = _family;
        int num_sons = l->num_sons();

        std::string info_file = _lib + "/RNA/records/" + (l->is_open() ? "open_" : "") + "loop";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        std::string line;
        int num = 0;
        while (std::getline(ifile, line) && set_loop_rec(templ_rec, line)) {
            if (std::find(m_disused_pdbs.begin(), m_disused_pdbs.end(), templ_rec._src) != m_disused_pdbs.end()) {
                continue;
            } else if (_is_virtual[l] && templ_rec._type == num_sons) {
                templ_rec._score = 0;
                if (templ_rec._src == _name) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if ((!(_source_pdb.empty())) && (templ_rec._src == _source_pdb)) {
                    _records[l].first.clear(); 
                    _records[l].first.push_back(templ_rec);
                    num = 1;
                    _find_self[l] = true;
                    break;
                }
            } else if (NucSS::pure_ss(NucSS::lower_ss(templ_rec._ss, 1)) == lower_ss) {
                templ_rec._score = (templ_rec._ss == ss ? 5 : 0);
                if (templ_rec._src == _name) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    if (_is_test) continue; else templ_rec._score += 2;
                }
                for (int i = 0; i < templ_rec._seq.size(); i++) {
                    if (seq[i] == templ_rec._seq[i]) {
                        if (ss[i] == templ_rec._ss[i] && ss[i] != '.' && ss[i] != '(' && ss[i] != ')') {
                            templ_rec._score += 2;
                        } else if (ss[i] == '(' || ss[i] == ')') {
                            templ_rec._score += 0.2;
                        } else {
                            templ_rec._score += 1;
                        }
                    }
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if ((!(_source_pdb.empty())) && (templ_rec._src == _source_pdb)) {
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
    }

    void find_helix_records(loop *l) {
        if (l->s.empty()) return;

        std::string seq = l->s.seq(), ss = l->s.ss(), family = _family;
        int len = l->s.len();

        std::string info_file = _lib + "/RNA/records/helix";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        std::string line;
        int num = 0;
        while (std::getline(ifile, line) && set_helix_rec(templ_rec, line)) {
            /* if (std::find(m_disused_pdbs.begin(), m_disused_pdbs.end(), templ_rec._src) != m_disused_pdbs.end()) {
                continue;
            } else */ if (templ_rec._len == len) {
                templ_rec._score = 0;
                if (templ_rec._src == jian::upper(_name)) {
                    templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    templ_rec._score += 2;
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
    }

};

} // namespace nuc3d
} // namespace jian

