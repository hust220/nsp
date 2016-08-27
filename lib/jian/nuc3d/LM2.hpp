#pragma once

#include "../pdb.hpp"
#include <nuc2d/N2D.h>
#include <dg/DG.h>
#include "../geom.hpp"
#include "Connect.hpp"
#include "Transform.hpp"
#include "JobPredict3D.hpp"
#include "../utils/log.hpp"

namespace jian {

class LM2 : public virtual JobPredict3D {
public:
    DG _dg;
    std::map<int, int> _m1, _m2;
    std::string _job;
    MatrixXf _native_helices;
    MatrixXf _helices;
    MatrixXf _native_scaffold;
    MatrixXf _scaffold;
    Model _native_model;

    std::map<string, MatrixXf> _mono_nuc_pars;
    std::map<string, MatrixXf> _adj_nuc_pars;
    std::map<string, MatrixXf> _aa_pars;

    N2D _ss_tree;
    std::vector<std::vector<res>> _helix_anchors;
    std::map<loop *, MatrixXf> _helix_coords;
    std::map<loop *, MatrixXf> _nuc_coords;
    std::map<loop *, MatrixXf> _atom_coords;
    MatrixXf _bound;
    DG dg;
    int _atom_nums_per_nuc = 5;
    double _err_radius = 0.001;
    Model _helix_file;
    int _view = 0;
    std::string _frag_par_path;
    std::vector<std::string> _frag_3_names;
    std::vector<std::string> _frag_4_names;
    std::vector<std::string> _frag_5_names;
    std::vector<std::array<double,  3>> _frag_3_dists;
    std::vector<std::array<double,  6>> _frag_4_dists;
    std::vector<std::array<double, 10>> _frag_5_dists;

    std::string _constraint_file;

//    LM2(const Par &);
//    Model operator ()();
//    void read_pars();
//    void set_native();
//    void set_bound();
//    void init();
//    Model run();
//    std::string parse_molecule_type(std::string seq);

//    void set_helix_anchors();
//    void set_helix_anchors(loop *src);
//    void set_helix_coords(loop *src);
//    Model to_all_atom(const MatrixXf &);
//    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> get_fragments();
//    MatrixXf make_helix_strand(MatrixXf, MatrixXf, int);
//    std::pair<MatrixXf, std::deque<Residue>> make_helix(const MatrixXf &, int);
//    std::pair<MatrixXf, std::deque<Residue>> make_frag(const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> &, MatrixXf, MatrixXf);

//    MatrixXf get_helix_par(const R5P &);
//    R5P get_helix(const helix &);
//    R5P create_helix(const string &);

//    void set_constraints();
//    MatrixXf apply_constraints(const MatrixXf &);
//    std::vector<std::vector<std::vector<int>>> read_constraints_file();
//    MatrixXf read_constraints_pos(std::string, int);
//    MatrixXf read_constraints_par(std::string, int);

    LM2(const Par &par) : JobPredict3D(par) {
        init();
    }

    Model operator ()() {
        return Transform(run())(_type, _seq);
    }

    // # Initialization
    void init() {
        read_pars();
        set_helix_anchors();
        set_native();
        set_bound();
    }

    // # Read parameters
    void read_pars() {
        // ## Read helix parameters
        std::string helix_file_name = _lib + "/RNA/pars/nuc3d/helix_coord.pdb";
        _helix_file = Model(helix_file_name);

        // ## Set the name of fragment parameter files
        _frag_par_path = _lib + "/RNA/pars/nuc3d/frag";
        std::string frag_3_par = _frag_par_path + "/frag-3.par";
        std::string frag_4_par = _frag_par_path + "/frag-4.par";
        std::string frag_5_par = _frag_par_path + "/frag-5.par";

        int num_frags = 0;
        // ## Read the information of fragments of 3 nt
        std::ifstream ifile(frag_3_par.c_str());
        ifile >> num_frags;
        _frag_3_names.resize(num_frags);
        _frag_3_dists.resize(num_frags);
        for (int i = 0; i < num_frags; i++) {
            ifile >> _frag_3_names[i];
            for (int j = 0; j < 3; j++) {
                ifile >> _frag_3_dists[i][j];
            }
        }
        ifile.close();

        // ## Read the information of fragments of 4 nt
        ifile.open(frag_4_par.c_str());
        ifile >> num_frags;
        _frag_4_names.resize(num_frags);
        _frag_4_dists.resize(num_frags);
        for (int i = 0; i < num_frags; i++) {
            ifile >> _frag_4_names[i];
            for (int j = 0; j < 6; j++) {
                ifile >> _frag_4_dists[i][j];
            }
        }
        ifile.close();

        // ## Read the information of fragments of 5 nt
        ifile.open(frag_5_par.c_str());
        ifile >> num_frags;
        _frag_5_names.resize(num_frags);
        _frag_5_dists.resize(num_frags);
        for (int i = 0; i < num_frags; i++) {
            ifile >> _frag_5_names[i];
            for (int j = 0; j < 10; j++) {
                ifile >> _frag_5_dists[i][j];
            }
        }
        ifile.close();
    }

    template<typename List> std::deque<Residue> find_best_frag_model(List &&list) {
        typedef struct {
            bool operator ()(const std::pair<int, double> &p1, const std::pair<int, double> &p2) {return p1.second < p2.second;}
        } Compare;
        std::priority_queue<std::pair<int, double>, std::deque<std::pair<int, double>>, Compare> rank;
        int size = std::distance(std::begin(list), std::end(list));
        if (size == 3) {
            for (int i = 0; i < _frag_3_dists.size(); i++) 
                rank.push(std::make_pair(i, geometry::distance(_frag_3_dists[i], std::forward<List>(list))));
            LOG << _frag_par_path + "/" + _frag_3_names[rank.top().first] + ".pdb" << std::endl;
            return get_residues_from_file(_frag_par_path + "/" + _frag_3_names[rank.top().first] + ".pdb");
        } else if (size == 10) {
            for (int i = 0; i < _frag_5_dists.size(); i++) 
                rank.push(std::make_pair(i, geometry::norm(_frag_5_dists[i], std::forward<List>(list), 10)));
            LOG << _frag_par_path + "/" + _frag_5_names[rank.top().first] + ".pdb" << std::endl;
            return get_residues_from_file(_frag_par_path + "/" + _frag_5_names[rank.top().first] + ".pdb");
        } else {
            throw "JIAN::NUC3D::LM2::find_best_frag_model error!";
        }
    }

    std::string parse_molecule_type(std::string seq) {
        std::regex nuc_pattern("[AUTGC]+");
        if (std::regex_match(seq, nuc_pattern)) {
            if (std::regex_search(seq, std::regex("T+"))) {
                return "DNA";
            } else {
                return "RNA";
            }
        } else {
            throw "JIAN::NUC3D::LM2::parse_molecule_type(std::string) error! Only support RNA and DNA.";
        }
    }

    void set_helix_anchors() {
        std::map<char, char> temp_map{{'(', '.'}, {')', '.'}, {'[', '.'}, {']', '.'}, 
                                      {'{', '.'}, {'}', '.'}, {'<', '.'}, {'>', '.'}, 
                                      {'A', '.'}, {'a', '.'}, {'B', '.'}, {'b', '.'}, 
                                      {'C', '.'}, {'c', '.'}, {'D', '.'}, {'d', '.'}, 
                                      {'E', '.'}, {'e', '.'}, {'F', '.'}, {'f', '.'}, 
                                      {'/', '.'}, {'.', '.'}, {'&', '&'}};
        for (auto &&pair : paired_keys) {
            if (std::count_if(_ss.begin(), _ss.end(), [&](const char &c){return c == pair.first || c == pair.second;})) {
                temp_map[pair.first] = '(';
                temp_map[pair.second] = ')';

                std::string str = _ss;
                std::transform(str.begin(), str.end(), str.begin(), [&](const char &c) {return temp_map[c];});
                N2D n2d2(str, _seq);
                set_helix_anchors(n2d2.head);

                temp_map[pair.first] = '.';
                temp_map[pair.second] = '.';;
            }
        }

        // ## Construct residue list and then sort
        std::vector<res> res_list;
        LOG << "helix anchors: " << std::endl;
        for (auto &&anchor: _helix_anchors) {
            for (auto &&r: anchor) {
                LOG << r.name << '-' << r.num << ' ';
                res_list.push_back(r);
            }
            LOG << std::endl;
        }
        std::sort(res_list.begin(), res_list.end(), [](const res &res1, const res &res2){return res1.num < res2.num;});

        // ## Construct map for the convenience of index convertion
        for (int i = 0; i < res_list.size(); i++) {
            _m1[res_list[i].num - 1] = i;    
            _m2[i] = res_list[i].num - 1;
        }

    }

    void set_helix_anchors(loop *src) {
        if (src == NULL) {
            return;
        } else {
            if (src->s.head != NULL) {
                std::vector<res> vec;
                bp *b = src->s.head;
                vec.push_back(b->res1);
                vec.push_back(b->res2);
                for (; b->next != NULL; b = b->next);
                vec.push_back(b->res1);
                vec.push_back(b->res2);
                _helix_anchors.push_back(vec);
            }
            set_helix_anchors(src->son);
            set_helix_anchors(src->brother);
        }
    }

    void set_native() {
        if (_native == "") {
            return;
        } else {
            _native_model = Model(_native);

            int len = _native_model.res_nums();
            if (len != _seq.size()) {
                throw "JIAN::NUC3D::LM2::set_native() error!";
            } else {
                _native_scaffold.resize(len, 3);
                int index = 0;
                for (auto &&chain: _native_model.chains) {
                    for (auto &&res: chain.residues) {
                        auto atom = res["C4*"].pos();
                        for (int j = 0; j < 3; j++) {
                            _native_scaffold(index, j) = atom[j];
                        }
                        index++;
                    }
                }
            }

            if (_helix_anchors.empty()) {
                return;
            } else {
                int len2 = _helix_anchors.size() * 4;
                _native_helices.resize(len2, 3);
                for (int i = 0; i < len2; i++) {
                    for (int j = 0; j < 3; j++) {
                        _native_helices(i, j) = _native_scaffold(_m2[i], j);
                    }
                }
            }
        }
    }

    void set_bound() {
        int len = _helix_anchors.size() * 4;
        if (len == 0) {
            return;    
        } else {
            _bound = MatrixXf::Zero(len, len);
            for (int i = 0; i < len; i++) {
                for (int j = i; j < len; j++) {
                    if (i == j) {
                        _bound(i, j) = 0;    
                    } else {
                        _bound(i, j) = 6.1 * (_m2[j] - _m2[i]);
                        _bound(j, i) = 6.1;
                    }
                }
            }

            // ## Add helix parameters into the bound matrix
            static auto fn_a = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 2 * 3.14159)) + (2.84 * n) * (2.84 * n));};
            static auto fn_c = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 1.5 * 3.14159)) + (2.84*n-4) * (2.84*n-4));};
            static auto fn_d = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 0.5 * 3.14159)) + (2.84*n+4) * (2.84*n+4));};
            for (auto &&anchor: _helix_anchors) {
                int i1 = _m1[anchor[0].num - 1];
                int i2 = _m1[anchor[1].num - 1];
                int i3 = _m1[anchor[2].num - 1];
                int i4 = _m1[anchor[3].num - 1];
                int len = _m2[i3] - _m2[i1];
                _bound(i1, i2) = _bound(i2, i1) = 15.1;
                _bound(i3, i4) = _bound(i4, i3) = 15.1;
                _bound(i1, i3) = _bound(i3, i1) = fn_a(len);
                _bound(i4, i2) = _bound(i2, i4) = fn_a(len);
                _bound(i1, i4) = _bound(i4, i1) = fn_c(len);
                _bound(i3, i2) = _bound(i2, i3) = fn_d(len);
            }

            // ## Calculate scaffold coordinates
            _dg = DG(_bound);
            LOG << "\nscaffold bound matrix: " << std::endl;
            LOG << _bound << "\n" << std::endl;
        }
    }

    Model run() {
        // ## Initialize variables
        int num_residues = _seq.size();
        _scaffold = MatrixXf::Zero(num_residues, 3);
        std::vector<Residue> residues(num_residues);

        if (!_helix_anchors.empty()) {

            _helices = _dg();
            LOG << "DG...\n";
            LOG << "scaffold energy: " << dg.E << std::endl;
            LOG << "scaffold: " << std::endl;
            LOG << _helices << std::endl;
            if (_native != "") LOG << "RMSD: " << geom::rmsd(_helices, _native_helices) << std::endl;

            // ## Copy scaffold coordinates to integral coordinates
            for (int i = 0; i < _helices.rows(); i++) {
                for (int j = 0; j < _helices.cols(); j++) {
                    _scaffold(_m2[i], j) = _helices(i, j);
                }
            }
            LOG << "\nIntegral coordinates: \n" << _scaffold << std::endl;

            // ## Calculate helix coordinates
            LOG << "\nCalculate helix coordinates:\n";
            for (auto &&anchor: _helix_anchors) {
                MatrixXf anchor_coord(4, 3); // the index order of an anchor is n0-n3 n1-n2
                for (int i = 0; i < 3; i++) {
                    anchor_coord(0, i) = _scaffold(anchor[0].num - 1, i);
                    anchor_coord(1, i) = _scaffold(anchor[2].num - 1, i);
                    anchor_coord(2, i) = _scaffold(anchor[3].num - 1, i);
                    anchor_coord(3, i) = _scaffold(anchor[1].num - 1, i);
                }
                int len = anchor[2].num - anchor[0].num + 1;

                // ### Construct helix
                auto helix_pair = make_helix(anchor_coord, len);
                auto &helix_coord = helix_pair.first;
                auto &helix_model = helix_pair.second;

                // ### Deposit helix information
                for (int i = 0; i < len; i++) {
                    int m = anchor[0].num - 1 + i;
                    int n = anchor[3].num - 1 + i;
                    for (int j = 0; j < 3; j++) {
                        _scaffold(m, j) = helix_coord(i, j);    
                        _scaffold(n, j) = helix_coord(i + len, j);
                    }
                    std::swap(residues[m], helix_model[i]);
                    std::swap(residues[n], helix_model[i + len]);
                }
                LOG << "superposed helix:\n" << helix_coord << std::endl;
            }
            LOG << "\nIntegral scaffold: \n" << _scaffold << std::endl;
        }

        // ## Construct fragments
        auto fragments = get_fragments();

        // ## Calculate coordinates of fragments
        MatrixXf coords;
        for (auto &&frag: fragments) {
            auto beg = std::get<0>(frag);
            auto center = std::get<1>(frag);
            auto end = std::get<2>(frag);
            MatrixXf a, b; // a is the matrix of beg, b is the matrix of end
            if (beg.size() != 0) {
                a.resize(2, 3);
                for (int i = 0; i < 3; i++) {
                    a(0, i) = _scaffold(beg[0], i);
                    a(1, i) = _scaffold(beg[1], i);
                }
                LOG << "\na: \n" << a << std::endl;
            }
            if (end.size() != 0) {
                b.resize(2, 3);
                for (int i = 0; i < 3; i++) {
                    b(0, i) = _scaffold(end[0], i);
                    b(1, i) = _scaffold(end[1], i);
                }
                LOG << "\nb: \n" << b << std::endl;
            }

            // ### Construct fragment
            auto frag_result = make_frag(frag, a, b);
            auto &frag_coord = frag_result.first;
            auto &frag_residues = frag_result.second;

            // ### Deposit the information of fragment
            int temp = 2;
            if (beg.size() == 0) temp = 0;
            for (int i = 0; i < center.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    _scaffold(center[i], j) = frag_coord(temp + i, j);
                }
                std::swap(residues[center[i]], frag_residues[i]);
            }
        }


        LOG << "\nfinal coords: \n" << _scaffold << std::endl;

        Model model;
        Chain chain;
        std::swap(chain.residues, residues);
        model.chains.push_back(chain);
        return model;

    }

    // # Make a helix fixed on two anchcors
    std::pair<MatrixXf, std::deque<Residue>> make_helix(const MatrixXf &anchor, int len) {
        int num_residues = len * 2;
        int num_residues_helix_file = _helix_file.res_nums();
        std::vector<int> sub_indices(num_residues);
        for (int i = 0; i < len; i++) {
            sub_indices[i] = i; 
            sub_indices[num_residues - 1 - i] = num_residues_helix_file - 1 - i;
        }
        auto residues = _helix_file.residues(sub_indices);

        MatrixXf scaffold(num_residues, 3);
        int num_residue = 0;
        for (int i = 0; i < residues.size(); i++) {
            auto pos = residues[i]["C4*"].pos();
            for (int j = 0; j < 3; j++)
                scaffold(i, j) = pos[j];
        }

        LOG << "helix:\n" << scaffold << std::endl;
        SupPos sp;
        sp(scaffold, hstack(scaffold.row(0), scaffold.row(len - 1), scaffold.row(len), scaffold.row(2 * len - 1)), anchor);
        auto c1 = -sp.c1;
        for (auto &&res: residues) {
            for (auto &&atom: res.atoms) {
                geom::move(atom, c1);
                geom::rotate(atom, sp.rot);
                geom::move(atom, sp.c2);
            }
        }
        return std::make_pair(scaffold, residues);
    }

    MatrixXf make_helix_strand(MatrixXf head, MatrixXf tail, int len) {
        MatrixXf bound(len + 1, len + 1);
        for (int i = 0; i < len + 1; i++) {
            for (int j = i; j < len + 1; j++) {
                if (i == j) {
                    bound(i, j) = 0;
                } else {
                    int d = j - i;
                    bound(i, j) = bound(j, i) = sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * d - 2 * 3.14159)) + (2.84*d) * (2.84*d));
                }
            }
        }
        DG dg(bound);
        auto coord = dg();
        return coord;
    }

    std::pair<MatrixXf, std::deque<Residue>> make_frag(const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> &frag, MatrixXf a, MatrixXf b) {
        auto &beg = std::get<0>(frag);
        auto &center = std::get<1>(frag);
        auto &end = std::get<2>(frag);
        std::vector<int> res_list;
        std::copy(std::get<0>(frag).begin(), std::get<0>(frag).end(), std::back_inserter(res_list));
        std::copy(std::get<1>(frag).begin(), std::get<1>(frag).end(), std::back_inserter(res_list));
        std::copy(std::get<2>(frag).begin(), std::get<2>(frag).end(), std::back_inserter(res_list));
        int len = res_list.size();
        std::map<int, int> m1, m2;
        for (int i = 0; i < len; i++) {
            m1[res_list[i]] = i;
            m2[i] = res_list[i];
        }

        MatrixXf bound = MatrixXf::Zero(len, len);
        for (int i = 0; i < len; i++) {
            for (int j = i; j < len; j++) {
                if (i == j) {
                    bound(i, j) = 0;    
                } else {
                    bound(i, j) = 6.1 * (m2[j] - m2[i]);
                    bound(j, i) = 6.1;
                }
            }
        }
        if (a.rows() == 0 && b.rows() == 0) {
            // pass
        } else if (a.rows() == 0) {
            double dist = (b.row(0) - b.row(1)).norm();
            bound(len - 2, len - 1) = dist;
            bound(len - 1, len - 2) = dist;
        } else if (b.rows() == 0) {
            double dist = (a.row(0) - a.row(1)).norm();
            bound(0, 1) = dist;
            bound(1, 0) = dist;
        } else {
            double dist1 = (a.row(0) - a.row(1)).norm();
            bound(0, 1) = bound(1, 0) = dist1;
            double dist2 = (a.row(1) - b.row(0)).norm();
            bound(1, len - 2) = bound(len - 2, 1) = dist2;
            double dist3 = (b.row(0) - b.row(1)).norm();
            bound(len - 2, len - 1) = bound(len - 1, len - 2) = dist3;
            double dist4 = (a.row(0) - b.row(1)).norm();
            bound(0, len - 1) = bound(len - 1, 0) = dist4;
            double dist5 = (a.row(0) - b.row(0)).norm();
            bound(0, len - 2) = bound(len - 2, 0) = dist5;
            double dist6 = (a.row(1) - b.row(1)).norm();
            bound(1, len - 1) = bound(len - 1, 1) = dist6;
        }
        LOG << "\nfragment bound: \n" << bound << std::endl;
        DG dg(bound);
        auto coord = dg();
    //    LOG << coord << "\n" << "energy: " << dg.E << "\n" << std::endl;
        LOG << "DG...\nenergy: " << dg.E << "\n" << coord << std::endl;
        LOG << "distances: " << std::endl;
        for (int i = 0; i < coord.rows() - 1; i++) {
            LOG << (coord.row(i) - coord.row(i + 1)).norm() << ' ';
        }
        LOG << std::endl;

        // ## Superpose matrix
        SupPos sp;
        if (a.rows() != 0 || b.rows() != 0) {
            int temp_len = a.rows() + b.rows();
            MatrixXf x(temp_len, 3), y(temp_len, 3);
            if (beg.size() != 0) {
                for (int i = 0; i < 3; i++) {
                    x(0, i) = coord(0, i);
                    x(1, i) = coord(1, i);
                    y(0, i) = a(0, i);
                    y(1, i) = a(1, i);
                }
            }
            if (end.size() != 0) {
                for (int i = 0; i < 3; i++) {
                    x(x.rows() - 2, i) = coord(coord.rows() - 2, i);
                    x(x.rows() - 1, i) = coord(coord.rows() - 1, i);
                    y(y.rows() - 2, i) = b(0, i);
                    y(y.rows() - 1, i) = b(1, i);
                }
            }
            sp(coord, x, y);
            LOG << "after superpose: " << std::endl;
            LOG << coord << std::endl;
        }

        // ## Find similar fragment model
        std::deque<Residue> residues;
        if (beg.empty() && end.empty()) {
            for (int i = 0; i < center.size() - 4; i++) {
                std::array<double, 10> dists;
                for (int j = 0, index = 0; j < 5; j++) {
                    for (int k = j + 1; k < 5; k++) {
                        dists[index] = (coord.row(i + j) - coord.row(i + k)).norm();
                        index++;
                    }
                }
                auto best_frag_model = find_best_frag_model(dists);
                MatrixXf x(5, 3), y(5, 3);
                for (int j = 0; j < 5; j++) {
                    for (int k = 0; k < 3; k++) {
                        x(j, k) = best_frag_model[j]["C4*"][k];
                        y(j, k) = coord(i + j, k);
                    }
                }
                sp(x, y);
                auto c1 = -sp.c1;
                for (auto &&res: best_frag_model) {
                        for (auto &&atom: res.atoms) {
                            geom::move(atom, c1);    
                            geom::rotate(atom, sp.rot);    
                            geom::move(atom, sp.c2);    
                        }
                }
                if (i == 0) {
                    residues.push_back(best_frag_model[0]);
                    residues.push_back(best_frag_model[1]);
                } 
                residues.push_back(best_frag_model[2]);
                if (i == center.size() - 5) {
                    residues.push_back(best_frag_model[3]);
                    residues.push_back(best_frag_model[4]);
                }
            }
        } else if (beg.empty()) {
            for (int i = center.size() - 1; i >= 0; i--) {
                auto best_frag_model = find_best_frag_model(std::vector<double>{(coord.row(i) - coord.row(i + 1)).norm(),
                                                             (coord.row(i) - coord.row(i + 2)).norm(),
                                                             (coord.row(i + 1) - coord.row(i + 2)).norm()});
                MatrixXf x(3, 3), y(3, 3);
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        x(j, k) = best_frag_model[j]["C4*"][k];
                        y(j, k) = coord(i + j, k);
                    }
                }
                sp(x, y);
                auto c1 = -sp.c1;
                for (auto &&res: best_frag_model) {
                        for (auto &&atom: res.atoms) {
                            geom::move(atom, c1);    
                            geom::rotate(atom, sp.rot);    
                            geom::move(atom, sp.c2);    
                        }
                }
                residues.push_back(best_frag_model[0]);
            }
        } else if (end.empty()) {
            for (int i = 0; i < center.size(); i++) {
                auto best_frag_model = find_best_frag_model(std::vector<double>{(coord.row(i) - coord.row(i + 1)).norm(),
                                                             (coord.row(i) - coord.row(i + 2)).norm(),
                                                             (coord.row(i + 1) - coord.row(i + 2)).norm()});
                MatrixXf x(3, 3), y(3, 3);
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        x(j, k) = best_frag_model[j]["C4*"][k];
                        y(j, k) = coord(i + j, k);
                    }
                }
                sp(x, y);
                auto c1 = -sp.c1;
                for (auto &&res: best_frag_model) {
                        for (auto &&atom: res.atoms) {
                            geom::move(atom, c1);    
                            geom::rotate(atom, sp.rot);    
                            geom::move(atom, sp.c2);    
                        }
                }
                residues.push_back(best_frag_model[2]);
            }
        } else {
            for (int i = 0; i < center.size(); i++) {
                std::array<double, 10> dists;
                for (int j = 0, index = 0; j < 5; j++) {
                    for (int k = j + 1; k < 5; k++) {
                        dists[index] = (coord.row(i + j) - coord.row(i + k)).norm();
                        index++;
                    }
                }
                auto best_frag_model = find_best_frag_model(dists);
                MatrixXf x(5, 3), y(5, 3);
                for (int j = 0; j < 5; j++) {
                    for (int k = 0; k < 3; k++) {
                        x(j, k) = best_frag_model[j]["C4*"][k];
                        y(j, k) = coord(i + j, k);
                    }
                }
                sp(x, y);
                auto c1 = -sp.c1;
                for (auto &&res: best_frag_model) {
                        for (auto &&atom: res.atoms) {
                            geom::move(atom, c1);    
                            geom::rotate(atom, sp.rot);    
                            geom::move(atom, sp.c2);    
                        }
                }
                residues.push_back(best_frag_model[2]);
            }
        }
        // ## Superpose fragment model
        return std::make_pair(coord, residues);
    }

    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> get_fragments() {
        std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> fragments;

        // ## Construct anchor list
        std::vector<int> anchor_list;
        for (auto &&vec: _helix_anchors) {
            for (auto &&res: vec) {
                anchor_list.push_back(res.num - 1);
            }
        }

        if (anchor_list.empty()) {
            std::vector<int> a, b(_seq.size()), c;
            std::iota(b.begin(), b.end(), 0);
            fragments.push_back(std::make_tuple(a, b, c));
        } else {
            std::sort(anchor_list.begin(), anchor_list.end(), [](int a, int b){return a < b;});

            // ## Construct the first fragment
            if (anchor_list[0] != 0) {
                std::vector<int> a, b, c; // a is the begin, b is the center, c is the end
                for (int j = 0; j < anchor_list[0]; j++) {
                    b.push_back(j);
                }
                c.push_back(anchor_list[0]);
                c.push_back(anchor_list[0] + 1);
                fragments.push_back(std::make_tuple(a, b, c));
            }
            // ## Construct the center fragments
            for (int i = 0; i < anchor_list.size() / 2 - 1; i++) {
                std::vector<int> a, b, c;
                a.push_back(anchor_list[i * 2 + 1] - 1);
                a.push_back(anchor_list[i * 2 + 1]);
                for (int j = anchor_list[i * 2 + 1] + 1; j < anchor_list[i * 2 + 2]; j++) {
                    b.push_back(j);
                }
                c.push_back(anchor_list[i * 2 + 2]);
                c.push_back(anchor_list[i * 2 + 2] + 1);
                fragments.push_back(std::make_tuple(a, b, c));
            }
            // ## Construct the last fragment
            if (anchor_list.back() != _seq.size() - 1) {
                std::vector<int> a, b, c;
                a.push_back(anchor_list.back() - 1);
                a.push_back(anchor_list.back());
                for (int j = anchor_list.back() + 1; j < _seq.size(); j++) {
                    b.push_back(j);
                }
                fragments.push_back(std::make_tuple(a, b, c));
            }

            // ## Print fragments information
            LOG << "\nfragments: " << std::endl;
            for (auto &&frag: fragments) {
                for (auto &&i: std::get<0>(frag)) LOG << i << '-';
                for (auto &&i: std::get<1>(frag)) LOG << i << '-';
                for (auto &&i: std::get<2>(frag)) LOG << i << '-';
                LOG << ' ';
                LOG << std::endl;
            }
        }
        return fragments;
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

    void set_helix_coords(loop *src) {
        for (bp *b = src->s.head; b != NULL; b = b->next) {
            if (std::count_if(_helix_anchors.begin(), _helix_anchors.end(), [&](const std::vector<res> &vec){
                return std::count_if(vec.begin(), vec.end(), [&](const res &r2){
                    return b->res1.type == r2.type && b->res1.name == r2.name && b->res1.num == r2.num;});})) {
                LOG << b->res1.num << ' ' << b->res1.type << ' ' << b->res1.name << std::endl;
            }
            if (std::count_if(_helix_anchors.begin(), _helix_anchors.end(), [&](const std::vector<res> &vec){return std::count_if(vec.begin(), vec.end(), [&](const res &r2){return b->res2.type == r2.type && b->res2.name == r2.name && b->res2.num == r2.num;});})) {
                LOG << b->res2.num << ' ' << b->res2.type << ' ' << b->res2.name << std::endl;
            }
        }
        for (res *r = src->head; r != NULL; r = r->next) {
            if (std::count_if(_helix_anchors.begin(), _helix_anchors.end(), [&](const std::vector<res> &vec){return std::count_if(vec.begin(), vec.end(), [&](const res &r2){return r->type == r2.type && r->name == r2.name && r->num == r2.num;});})) {
                LOG << r->num << ' ' << r->type << ' ' << r->name << std::endl;
            }
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

