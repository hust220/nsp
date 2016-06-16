/**
 * @file LoopModelling22.cpp
 * @brief Modelling loop by using distance geometry algorithm.
 *
 * Each nucleotide is represented by 5 atoms.
 *   A: C5* O3* C1* N6 C2
 *   U: C5* O3* C1* O2 O4
 *   G: C5* O3* C1* O6 N2
 *   C: C5* O3* C1* O2 N4
 *
 * @author wj_hust08@163.com
 * @version 1.0
 * @date 2015-5-25
*/

#pragma once

#include "../pdb.hpp"
#include <nuc2d/N2D.h>
#include <dg/DG.h>
#include "Connect.hpp"

namespace jian {

class LoopModelling2 {
public:
    using namespace nuc2d;

    std::string _seq;
    std::string _ss;
    MatrixXf _bound;
    DG dg;
    int _atom_nums_per_nuc = 5;
    double _err_radius = 0.001;
    std::string _lib;
    std::string _type = "RNA";
    int _view = 0;

    std::string _constraint_file;

//    LoopModelling2(string type = "RNA");
//    Model operator ()(string seq, string ss, string constraint_file = "", int num = 1);
//    Model to_all_atom(const MatrixXf &);

//    std::map<string, MatrixXf> _mono_nuc_pars;
//    std::map<string, MatrixXf> _adj_nuc_pars;
//    std::map<string, MatrixXf> _aa_pars;

//    void init();
//    void set_base_pairs(loop *);

//    MatrixXf get_helix_par(const R5P &);
//    R5P get_helix(const helix &);
//    R5P create_helix(const string &);

//    void set_constraints();
//    MatrixXf apply_constraints(const MatrixXf &);
//    std::vector<std::vector<std::vector<int>>> read_constraints_file();
//    MatrixXf read_constraints_pos(std::string, int);
//    MatrixXf read_constraints_par(std::string, int);

    LoopModelling2(std::string type = "RNA"): _type(type) {
        char *lib = getenv("NSP");
        assert(lib);
        _lib = string() + lib;
        _lib += "/" + boost::to_upper_copy(_type);

        /// read mononucleotide parameters
        string file_name = _lib + "/pars/5p/mono_nuc_pars";
        ifstream ifile(file_name.c_str());
        assert(ifile);
        for (int ii = 0; ii < 4; ii++) {
            string c;
            ifile >> c;
            _mono_nuc_pars[c].resize(_atom_nums_per_nuc, _atom_nums_per_nuc);
            for (int i = 0; i < _atom_nums_per_nuc; i++) {
                for (int j = 0; j < _atom_nums_per_nuc; j++) {
                    ifile >> _mono_nuc_pars[c](i, j);
                }
            }
        }
        ifile.close();

        /// read adjacent nucleotides parameters
        file_name = _lib + "/pars/5p/adj_nuc_pars";
        ifile.open(file_name.c_str());
        for (int ii = 0; ii < 16; ii++) {
            string str;
            ifile >> str;
            _adj_nuc_pars[str].resize(_atom_nums_per_nuc, _atom_nums_per_nuc);
            for (int i = 0; i < _atom_nums_per_nuc; i++) {
                for (int j = 0; j < _atom_nums_per_nuc; j++) {
                    ifile >> _adj_nuc_pars[str](i, j);
                }
            }
        }
        ifile.close();
    }

    Model operator ()(string seq, string ss, string constraint_file, int num) {
        _seq = seq;
        _ss = ss;
        _constraint_file = constraint_file;

        assert(seq.size() == count_if(ss.begin(), ss.end(), [](char c) {
            return set<char>{'.', '(', ')', '[', ']'}.count(c);
        }));

        init();
        dg = DG(_bound, 2);
        dg.log.set_display(_view);
        return to_all_atom(apply_constraints(dg()));
    }

    Model to_all_atom(const MatrixXf &scaffold) {
        int res_nums = _seq.size();
        Chain chain;

        BuildNuc build_nuc(_type);
        build_nuc._view = _view;
        for (int i = 0; i < res_nums; i++) {
            chain.residues.push_back(
                build_nuc(
                    (_type == "DNA" ? (string("D") + _seq[i]) : (string() + _seq[i])),
                    scaffold.block(i * _atom_nums_per_nuc, 0, _atom_nums_per_nuc, 3)
                )
            );
        }

        Model model;
        model.chains.push_back(chain);
        return model;
    }

    void init() {
        int res_nums = _seq.size();
        int atom_nums = res_nums * _atom_nums_per_nuc;

        /// begin initializing bound matrix
        _bound.resize(atom_nums, atom_nums);
        for (int i = 0; i < atom_nums; i++) {
            for (int j = i; j < atom_nums; j++) {
                if (i == j) {
                    _bound(i, j) = 0;
                } else {
                    _bound(j, i) = 2;
                    _bound(i, j) = 999;
                }
            }
        }
        /// end initializing bound matrix

        /// begin mononucleotide
        for (int i = 0; i < res_nums; i++) {
            for (int j = 0; j < _atom_nums_per_nuc; j++) {
                for (int k = j + 1; k < _atom_nums_per_nuc; k++) {
                    _bound(i * _atom_nums_per_nuc + j, i * _atom_nums_per_nuc + k)
                        = _mono_nuc_pars[string() + _seq[i]](j, k);
                    _bound(i * _atom_nums_per_nuc + k, i * _atom_nums_per_nuc + j)
                        = _mono_nuc_pars[string() + _seq[i]](j, k);
                }
            }
        }
        /// end mononucleotide

        /// begin adjacent nucleotides
        for (int i = 0; i + 1 < res_nums; i++) {
            int j = i + 1;
            for (int m = 0; m < _atom_nums_per_nuc; m++) {
                for (int n = 0; n < _atom_nums_per_nuc; n++) {
                    if (m != 1 || n != 0) continue;
                    string str = string() + _seq[i] + _seq[j];
                    _bound(i * _atom_nums_per_nuc + m, j * _atom_nums_per_nuc + n)
                        = _adj_nuc_pars[str](m, n);
                    _bound(j * _atom_nums_per_nuc + n, i * _atom_nums_per_nuc + m)
                        = _adj_nuc_pars[str](m, n);
                }
            }
        }
        /// end adjacent nucleotides

        /// begin base pairs
        N2D mol2d(_ss, _seq);
        set_base_pairs(mol2d.head);

        string str = _ss;
        transform(str.begin(), str.end(), str.begin(), [](const char &c) {
            std::map<char, char> temp_map{{'(', '['}, {')', ']'}, {'[', '('}, 
                                          {']', ')'}, {'.', '.'}, {'&', '&'}};
            return temp_map[c];
        });
        N2D mol2d_2(str, _seq);
        set_base_pairs(mol2d_2.head);
        /// end base pairs

        /// begin constraints
        set_constraints();
        /// end constraints

        /// chirality
    //    int flag = 0;
    //    for (int i = 3; i < len;) {
    //        if (ss[i / 6] == '(' && i % 6 == 5 && i != 5) {
    //            i += 4;
    //        } else {
    //            i++;
    //        }
    //        flag++;
    //    }
    //    delete chir;
    //    chir = new Matr_(flag, 5);
    //    for (int i = 3, flag = 0; i < len;) {
    //        chir->data[flag][0] = i - 3;
    //        chir->data[flag][1] = i - 2;
    //        chir->data[flag][2] = i - 1;
    //        chir->data[flag][3] = i;
    //        chir->data[flag][4] = chirality[(i - 2) % 6];
    //        if (ss[i / 6] == '(' && i % 6 == 5 && i != 5) {
    //            i += 4;
    //        } else {
    //            i++;
    //        }
    //        flag++;
    //    }

    }

    void set_base_pairs(loop *src) {
        if (src == NULL) {
            return;
        } else {
            if (src->s.head != NULL) {
                MatrixXf helix_par = get_helix_par(get_helix(src->s));
                int helix_len = src->s.getLen();
                int orig_pos_chain1 = (src->s.head->res1.num - 1) * _atom_nums_per_nuc;
                int orig_pos_chain2 = (src->s.head->res2.num - helix_len) * _atom_nums_per_nuc;
                int atom_nums_per_chain = _atom_nums_per_nuc * helix_len;
                for (int i: {0, 1}) {
                    for (int j: {0, 1}) {
                        for (int ii = 0; ii < atom_nums_per_chain; ii++) {
                            for (int jj = 0; jj < atom_nums_per_chain; jj++) {
                                if (i == j && ii == jj) {
                                    continue;
                                }
                                int a = orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * i + ii;
                                int b = orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * j + jj;
                                int c = i * (helix_par.rows() - atom_nums_per_chain) + ii;
                                int d = j * (helix_par.rows() - atom_nums_per_chain) + jj;
                                double dist = helix_par(c, d);
                                if (a < b) {
                                    _bound(a, b) = dist + _err_radius;
                                } else {
                                    _bound(a, b) = dist - _err_radius;
                                }
                            }
                        }
                    }
                }
            }
            set_base_pairs(src->son);
            set_base_pairs(src->brother);
        }
    }

    Model get_helix(const helix &h) {
        string file_name = _lib + "/info/helix";
        ifstream ifile(file_name.c_str());
        assert(ifile);
        string helix_name, helix_seq, helix_ss, helix_src;
        int helix_len;
        string target_seq = h.seq();
        string pdb_name;
        while (ifile) {
            ifile >> helix_name >> helix_len >> helix_seq >> helix_ss >> helix_src;
            if (target_seq == helix_seq) {
                pdb_name = _lib + "/helix/" + helix_name + ".pdb";
                break;
            }
        }
        ifile.close();
        if (pdb_name == "") {
            return create_helix(target_seq);
        } else {
            return R5P(pdb_name);
        }
    }

    Model create_helix(const string &seq) {
        if (seq.size() < 4) {
            std::cerr << "Assemble::createHelix error! The length of the helix"
                      << " to be create should not be less than 4!" << std::endl;
            exit(1);
        } else if (seq.size() % 2 == 1) {
            std::cerr << "Assemble::createHelix error! The length of the helix"
                      << " to be create should be even!" << std::endl;
            exit(1);
        } else if (seq.size() == 4 || seq.size() == 6) {
            string file_name = _lib + "/basepair/" + seq + ".pdb";
            ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XXXXXX.pdb";
            }
            ifile.close();
            return Model(file_name);
        } else {
            string file_name = _lib + "/basepair/" + seq.substr(0, 3) + 
                               seq.substr(seq.size() - 3, 3) + ".pdb";
            ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XXXXXX.pdb";
            }
            ifile.close();
            return Connect()(Model(file_name), create_helix(seq.substr(1, seq.size() - 2)), 2, 3);
        }
    }

    MatrixXf get_helix_par(const Model &r5p) {
        int atom_nums = r5p.atom_nums();
        MatrixXf helix_par(atom_nums, atom_nums);
        int num_1 = 0;
        for (auto &chain: r5p.chains) {
            for (auto &residue: chain.residues) {
                for (auto &atom: residue.atoms) {
                    int num_2 = 0;
                    for (auto &chain2: r5p.chains) {
                        for (auto &residue2: chain2.residues) {
                            for (auto &atom2: residue2.atoms) {
                                helix_par(num_1, num_2) = Point(atom.x, atom.y, atom.z).dist(
                                        Point(atom2.x, atom2.y, atom2.z));
                                num_2++;
                            }
                        }
                    }
                    num_1++;
                }
            }
        }
        return helix_par;
    }

    void set_constraints() {
        if (_constraint_file.empty()) return;

        /// Read constraints
        auto constraints = read_constraints_file();

        /// Set _bound
        for (auto &&temp_helix: constraints) {
            int layer_size = temp_helix[0].size();
            int layer_nums = temp_helix.size();
            auto par = read_constraints_par(std::map<int, std::string>{
                {4, "quadruplex"}, {3, "triple"}, {2, "pair"}
            }[layer_size], layer_nums);

            int index = 0;
            std::map<int, int> temp_map;
            for (auto &&vec: temp_helix) {
                for (int i = 0; i < layer_size * 5; i++) {
                    temp_map[index] = vec[i / 5] * 5 + i % 5;
                    index++;
                }
            }

            for (int i = 0; i < index; i++) {
                for (int j = 0; j < index; j++) {
    //                if (i % 5 > 1 && j % 5 > 1) {
                    if (i / 5 != j / 5) {
                        _bound(temp_map[i], temp_map[j]) = par(i, j);
                    }
                }
            }
        }
    }

    MatrixXf apply_constraints(const MatrixXf &scaffold) {
        if (_constraint_file.empty()) return scaffold;

        /// Read constraints
        auto constraints = read_constraints_file();

        /// Set _bound
        MatrixXf coords = scaffold;
        for (auto &&temp_helix: constraints) {
            int layer_size = temp_helix[0].size();
            int layer_nums = temp_helix.size();
            auto pos = read_constraints_pos(std::map<int, std::string>{
                {4, "quadruplex"}, {3, "triple"}, {2, "pair"}
            }[layer_size], layer_nums);

            int index = 0;
            std::map<int, int> temp_map;
            for (auto &&vec: temp_helix) {
                for (int i = 0; i < layer_size * 5; i++) {
                    temp_map[index] = vec[i / 5] * 5 + i % 5;
                    index++;
                }
            }

            MatrixXf source_mat = pos;
            MatrixXf target_mat(index, 3);

            for (int i = 0; i < index; i++) {
                for (int j = 0; j < 3; j++) {
                    target_mat(i, j) = coords(temp_map[i], j);
                }
            }

            SupPos()(pos, source_mat, target_mat);

            for (int i = 0; i < index; i++) {
                for (int j = 0; j < 3; j++) {
                    coords(temp_map[i], j) = pos(i, j);
                }
            }
        }

        return coords;
    }

    std::vector<std::vector<std::vector<int>>> read_constraints_file() {
        std::string constraint_type;
        std::vector<std::vector<std::vector<int>>> constraints;
        std::ifstream ifile(_constraint_file.c_str());
        std::string line;
        std::vector<std::vector<int>> temp_helix;
        int constraint_size;
        while (true) {
            if (!std::getline(ifile, line)) {
                if (!temp_helix.empty()) {
                    constraints.push_back(temp_helix);
                    temp_helix.clear();
                }
                break;
            }
            boost::trim(line);
            std::vector<std::string> splited_line;
            boost::split(splited_line, line, boost::is_any_of(" ,-"), boost::token_compress_on);
            if (splited_line == std::vector<std::string>{"G", "quadruplex", "helix"}) {
                constraint_size = 4;
                if (!temp_helix.empty()) {
                    constraints.push_back(temp_helix);
                    temp_helix.clear();
                }
                continue;
            } else if (splited_line == std::vector<std::string>{"base", "triple", "helix"}) {
                constraint_size = 3;
                if (!temp_helix.empty()) {
                    constraints.push_back(temp_helix);
                    temp_helix.clear();
                }
                continue;
            } else if (splited_line == std::vector<std::string>{"base", "pair", "helix"}) {
                constraint_size = 2;
                if (!temp_helix.empty()) {
                    constraints.push_back(temp_helix);
                    temp_helix.clear();
                }
                continue;
            }
            if (splited_line.size() != constraint_size) continue;
            std::vector<int> temp_vec;
            std::transform(std::begin(splited_line), std::end(splited_line), 
                           std::back_inserter(temp_vec), [](const std::string &s){
                return stoi(s) - 1;
            });
            temp_helix.push_back(temp_vec);
        }
        ifile.close();
        
        return constraints;
    }

    MatrixXf read_constraints_pos(std::string type, int layer_nums) {
        int layer_size = std::map<std::string, int>{{"quadruplex", 4}, {"triple", 3}, {"pair", 2}}[type];
        int atom_nums = 5 * layer_nums * layer_size;

        /// Read coordinates
        MatrixXf coords(atom_nums, 3);
        string par_file = _lib + "/pars/5p/" + type + ".pos";
        ifstream ifile(par_file.c_str());
        for (int i = 0; i < atom_nums; i++) {
            for (int j = 0; j < 3; j++) {
                ifile >> coords(i, j);
            }
        }
        ifile.close();

        return coords;
    }

    MatrixXf read_constraints_par(std::string type, int layer_nums) {
        auto coords = read_constraints_pos(type, layer_nums);
        int atom_nums = coords.rows();
        MatrixXf par(atom_nums, atom_nums);
        for (int i = 0; i < atom_nums; i++) {
            for (int j = 0; j < atom_nums; j++) {
                par(i, j) = (coords.row(i) - coords.row(j)).norm();
            }
        }

        return par;
    }

};

} // namespace jian

