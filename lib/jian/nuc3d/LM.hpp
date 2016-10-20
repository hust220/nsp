/**
 * @file LM.cpp
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
 * @date 2015-10-16
*/

#pragma once

#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "../dg/DG.hpp"
#include "Connect.hpp"
#include "BuildNuc.hpp"

namespace jian {

class LM {
public:
    using Mat = MatrixXd;

    std::map<string, Mat> _mono_nuc_pars;
    std::map<string, Mat> _adj_nuc_pars;
    std::map<string, Mat> _aa_pars;

    std::string _seq;
    std::string _ss;
    int _num_atoms_per_nuc = 5;
    double _err_radius = 0.001;
    std::string _lib;
    std::string _type = "RNA";
    int _view = 0;

    std::string _constraint_file;

    LM(std::string type = "RNA"): _type(type) {
        char *lib = getenv("NSP");
        assert(lib);
        _lib = string() + lib;
        _lib += "/" + jian::to_upper_copy(_type);

        /// read mononucleotide parameters
        std::string file_name = _lib + "/pars/5p/mono_nuc_pars";
        std::ifstream ifile(file_name.c_str());
        assert(ifile);
        for (int ii = 0; ii < 4; ii++) {
            string c;
            ifile >> c;
            _mono_nuc_pars[c].resize(_num_atoms_per_nuc, _num_atoms_per_nuc);
            for (int i = 0; i < _num_atoms_per_nuc; i++) {
                for (int j = 0; j < _num_atoms_per_nuc; j++) {
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
            _adj_nuc_pars[str].resize(_num_atoms_per_nuc, _num_atoms_per_nuc);
            for (int i = 0; i < _num_atoms_per_nuc; i++) {
                for (int j = 0; j < _num_atoms_per_nuc; j++) {
                    ifile >> _adj_nuc_pars[str](i, j);
                }
            }
        }
        ifile.close();
    }

    std::vector<Model> operator ()(string seq, string ss, std::string constraint_file = "", int num = 1) {
        _seq = seq; _ss = ss; _constraint_file = constraint_file;
        if (NucSS::len_ss(ss) != seq.size()) throw "jian::LM error!";
        DG dg(init()); dg.log.set_display(_view);
        std::vector<Model> models;
        for (int i = 0; i < num; i++) models.push_back(to_all_atom(apply_constraints(dg())));
        return models;
    }

    Model to_all_atom(const Mat &scaffold) {
        int num_residues = _seq.size();
        Chain chain;

        BuildNuc build_nuc(_type);
        build_nuc._view = _view;
        for (int i = 0; i < num_residues; i++) {
            chain.push_back(
                build_nuc(
                    (_type == "DNA" ? (string("D") + _seq[i]) : (string() + _seq[i])),
                    scaffold.block(i * _num_atoms_per_nuc, 0, _num_atoms_per_nuc, 3)
                )
            );
        }

        Model model;
        model.push_back(chain);
        return model;
    }

    Mat init() {
        int num_residues = _seq.size();
        int num_atoms = num_residues * _num_atoms_per_nuc;

        // Initialize bound matrix
        Mat bound = Mat::Zero(num_atoms, num_atoms);
        for (int i = 0; i < num_atoms; i++) for (int j = i + 1; j < num_atoms; j++) {
            if (i != j) { bound(j, i) = 2; bound(i, j) = 999; }
        }

        // Apply parameters of mono-nucleotide
        for (int i = 0; i < num_residues; i++) for (int j = 0; j < _num_atoms_per_nuc; j++) for (int k = j + 1; k < _num_atoms_per_nuc; k++) {
            bound(i * _num_atoms_per_nuc + j, i * _num_atoms_per_nuc + k) = _mono_nuc_pars[string() + _seq[i]](j, k);
            bound(i * _num_atoms_per_nuc + k, i * _num_atoms_per_nuc + j) = _mono_nuc_pars[string() + _seq[i]](j, k);
        }

        // Apply parameters of adjacent nucleotides
        for (int i = 0; i + 1 < num_residues; i++) {
            int j = i + 1;
            for (int m = 0; m < _num_atoms_per_nuc; m++) {
                for (int n = 0; n < _num_atoms_per_nuc; n++) {
                    if (m != 1 || n != 0) continue;
                    std::string str = std::string() + _seq[i] + _seq[j];
                    bound(i * _num_atoms_per_nuc + m, j * _num_atoms_per_nuc + n) = _adj_nuc_pars[str](m, n);
                    bound(j * _num_atoms_per_nuc + n, i * _num_atoms_per_nuc + m) = _adj_nuc_pars[str](m, n);
                }
            }
        }

        // Apply parameters of base pairs
        SSTree ss_tree1; ss_tree1.make(_seq, _ss); set_base_pairs(ss_tree1.head, bound);
        std::string str = _ss;
        transform(str.begin(), str.end(), str.begin(), [](const char &c) {
            std::map<char, char> temp_map{{'(', '['}, {')', ']'}, {'[', '('}, {']', ')'}, {'.', '.'}, {'&', '&'}};
            return temp_map[c];
        });
        SSTree ss_tree2; ss_tree2.make(_seq, str); set_base_pairs(ss_tree2.head, bound);

        set_constraints(bound);

        return bound;
    }

    void set_base_pairs(loop *src, Mat &bound) {
        if (src == NULL) {
            return;
        } else {
            if (src->s.head != NULL) {
                Mat helix_par = get_helix_par(get_helix(src->s));
                int helix_len = src->s.len();
                int orig_pos_chain1 = (src->s.head->res1.num - 1) * _num_atoms_per_nuc;
                int orig_pos_chain2 = (src->s.head->res2.num - helix_len) * _num_atoms_per_nuc;
                int atom_nums_per_chain = _num_atoms_per_nuc * helix_len;
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
                                    bound(a, b) = dist + _err_radius;
                                } else {
                                    bound(a, b) = dist - _err_radius;
                                }
                            }
                        }
                    }
                }
            }
            set_base_pairs(src->son, bound);
            set_base_pairs(src->brother, bound);
        }
    }

    Model get_helix(const helix &h) {
        std::string file_name = _lib + "/records/helix";
        std::ifstream ifile(file_name.c_str());
        assert(ifile);
        std::string helix_name, helix_seq, helix_ss, helix_src;
        int helix_len;
        std::string target_seq = h.seq();
        std::string pdb_name;
        while (ifile) {
            ifile >> helix_name >> helix_len >> helix_seq >> helix_ss >> helix_src;
            if (target_seq == helix_seq) {
                pdb_name = _lib + "/templates/" + helix_name + ".pdb";
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

    Model create_helix(const std::string &seq) {
        if (seq.size() < 4) {
            std::cerr << "Assemble::createHelix error! The length of the helix"
                      << " to be create should not be less than 4!" << std::endl;
            exit(1);
        } else if (seq.size() % 2 == 1) {
            std::cerr << "Assemble::createHelix error! The length of the helix"
                      << " to be create should be even!" << std::endl;
            exit(1);
        } else if (seq.size() == 4 || seq.size() == 6) {
            std::string file_name = _lib + "/basepair/" + seq + ".pdb";
            ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XXXXXX.pdb";
            }
            ifile.close();
            return Model(file_name);
        } else {
            std::string file_name = _lib + "/basepair/" + seq.substr(0, 3) + 
                               seq.substr(seq.size() - 3, 3) + ".pdb";
            ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XXXXXX.pdb";
            }
            ifile.close();
            return Connect()(Model(file_name), create_helix(seq.substr(1, seq.size() - 2)), 2, 3);
        }
    }

    Mat get_helix_par(const Model &r5p) {
        int atom_nums = num_atoms(r5p);
        Mat helix_par(atom_nums, atom_nums);
        int num_1 = 0;
        for (auto &chain: r5p) {
            for (auto &residue: chain) {
                for (auto &atom: residue) {
                    int num_2 = 0;
                    for (auto &chain2: r5p) {
                        for (auto &residue2: chain2) {
                            for (auto &atom2: residue2) {
                                helix_par(num_1, num_2) = Point(atom[0], atom[1], atom[2]).dist(
                                        Point(atom2[0], atom2[1], atom2[2]));
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

    void set_constraints(Mat &bound) {
        if (_constraint_file.empty()) return;

        auto constraints = read_constraints_file();
        std::map<int, std::string> types {{4, "quadruplex"}, {3, "triple"}, {2, "pair"}};
        for (auto &&temp_helix: constraints) {
            int layer_size = temp_helix[0].size(), layer_nums = temp_helix.size();
            auto par = read_constraints_par(types[layer_size], layer_nums);

            int index = 0; std::map<int, int> temp_map;
            for (auto &&vec: temp_helix) for (int i = 0; i < layer_size * 5; i++) {
                temp_map[index] = vec[i / 5] * 5 + i % 5;
                index++;
            }

            for (int i = 0; i < index; i++) for (int j = 0; j < index; j++) {
                if (i / 5 != j / 5) bound(temp_map[i], temp_map[j]) = par(i, j);
            }
        }
    }

    Mat apply_constraints(const Mat &scaffold) {
        if (_constraint_file.empty()) return scaffold;

        auto constraints = read_constraints_file();
        Mat coords = scaffold;
        std::map<int, std::string> types {{4, "quadruplex"}, {3, "triple"}, {2, "pair"}};
        for (auto &&temp_helix: constraints) {
            int layer_size = temp_helix[0].size(), layer_nums = temp_helix.size();
            auto pos = read_constraints_pos(types[layer_size], layer_nums);

            int index = 0; std::map<int, int> temp_map;
            for (auto &&vec: temp_helix) for (int i = 0; i < layer_size * 5; i++) {
                temp_map[index] = vec[i / 5] * 5 + i % 5;
                index++;
            }

            Mat source_mat = pos; Mat target_mat(index, 3);
            for (int i = 0; i < index; i++) for (int j = 0; j < 3; j++) target_mat(i, j) = coords(temp_map[i], j);
            geom::suppos(pos, source_mat, target_mat);
            for (int i = 0; i < index; i++) for (int j = 0; j < 3; j++) coords(temp_map[i], j) = pos(i, j);
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
            jian::trim(line);
            std::vector<std::string> splited_line;
            jian::tokenize(line, splited_line, " ,-");
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
                return std::stoi(s) - 1;
            });
            temp_helix.push_back(temp_vec);
        }
        ifile.close();
        
        return constraints;
    }

    Mat read_constraints_pos(std::string type, int layer_nums) {
        int layer_size = std::map<std::string, int>{{"quadruplex", 4}, {"triple", 3}, {"pair", 2}}[type];
        int atom_nums = 5 * layer_nums * layer_size;

        /// Read coordinates
        Mat coords(atom_nums, 3);
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

    Mat read_constraints_par(std::string type, int layer_nums) {
        auto coords = read_constraints_pos(type, layer_nums);
        int atom_nums = coords.rows();
        Mat par(atom_nums, atom_nums);
        for (int i = 0; i < atom_nums; i++) {
            for (int j = 0; j < atom_nums; j++) {
                par(i, j) = (coords.row(i) - coords.row(j)).norm();
            }
        }

        return par;
    }

};

} // namespace jian
