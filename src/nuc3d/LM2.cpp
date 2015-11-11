#include "LM2.h"
#include "BuildNuc.h"

namespace jian {
namespace nuc3d {

using namespace nuc2d;

LM2::LM2(string type): _type(type) {
}

std::vector<Model> LM2::operator ()(string seq, string ss, string constraint_file, int num) {
    _seq = seq;
    _ss = ss;
    _constraint_file = constraint_file;

    assert(seq.size() == count_if(ss.begin(), ss.end(), [](char c) {
        return set<char>{'.', '(', ')', '[', ']'}.count(c);
    }));

    init();
    std::vector<Model> models;
//    dg = DG(_bound, 2);
//    dg.view = _view;
//    std::vector<Model> models;
//    for (int i = 0; i < num; i++) {
//        models.push_back(to_all_atom(apply_constraints(dg())));
//    }
    return models;
}

Model LM2::to_all_atom(const MatrixXf &scaffold) {
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

void LM2::init() {
    /// Construct helix anchors
    _ss_tree = N2D(_ss, _seq);
    set_helix_anchors(_ss_tree.head);

    string str = _ss;
    map<char, char> temp_map{{'(', '['}, {')', ']'}, {'[', '('}, {']', ')'}, {'.', '.'}, {'&', '&'}};
    transform(str.begin(), str.end(), str.begin(), [&](const char &c) {return temp_map[c];});
    N2D n2d2(str, _seq);
    set_helix_anchors(n2d2.head);

    std::vector<res> res_list;
    std::cout << "helix anchors: " << std::endl;
    for (auto &&anchor: _helix_anchors) {
        for (auto &&r: anchor) {
            std::cout << r.name << '-' << r.num << ' ';
            res_list.push_back(r);
        }
        std::cout << std::endl;
    }
    std::sort(res_list.begin(), res_list.end(), [](const res &res1, const res &res2){return res1.num < res2.num;});

    /// Construct map
    std::map<int, int> m1, m2;
    for (int i = 0; i < res_list.size(); i++) {
        m1[res_list[i].num] = i;    
        m2[i] = res_list[i].num;
    }

    /// Initiate bound matrix
    int len = _helix_anchors.size() * 4;
    _bound = MatrixXf::Zero(len, len);
    _bound.resize(len, len);
    for (int i = 0; i < len; i++) {
        for (int j = i; j < len; j++) {
            if (i == j) {
                _bound(i, j) = 0;    
            } else {
                _bound(i, j) = 6.1 * (m2[j] - m2[i]);
                _bound(j, i) = 6.1;
            }
        }
    }

    /// Add helix parameters into the bound matrix
    auto fn_a = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 2 * 3.14159)) + (2.84 * n) * (2.84 * n));};
    auto fn_c = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 1.5 * 3.14159)) + (2.84*n-4) * (2.84*n-4));};
    auto fn_d = [](int n){return sqrt(2 * 9.7 * 9.7 * (1 - cos(0.562 * n - 0.5 * 3.14159)) + (2.84*n+4) * (2.84*n+4));};
    for (auto &&anchor: _helix_anchors) {
        int i1 = m1[anchor[0].num];
        int i2 = m1[anchor[1].num];
        int i3 = m1[anchor[2].num];
        int i4 = m1[anchor[3].num];
        int len = m2[i3] - m2[i1];
        _bound(i1, i2) = _bound(i2, i1) = 15.1;
        _bound(i3, i4) = _bound(i4, i3) = 15.1;
        _bound(i1, i3) = _bound(i3, i1) = fn_a(len);
        _bound(i4, i2) = _bound(i2, i4) = fn_a(len);
        _bound(i1, i4) = _bound(i4, i1) = fn_c(len);
        _bound(i3, i2) = _bound(i2, i3) = fn_d(len);
    }

    /// Calculate scaffold coordinates
    DG dg(_bound);
    std::cout << "\nscaffold bound matrix: " << std::endl;
    std::cout << _bound << "\n" << std::endl;
    _scaffold = dg();
    std::cout << "DG...\n";
    std::cout << "scaffold energy: " << dg.E << std::endl;
    std::cout << "scaffold: " << std::endl;
    std::cout << _scaffold << std::endl;

    /// Calculate helix coordinates
    std::cout << "\nCalculate helix coordinates:\n";
    std::map<int, MatrixXf> helix_coords;
    SupPos sp;
    for (auto &&anchor: _helix_anchors) {
        MatrixXf anchor_coord(4, 3);
        for (int i = 0; i < 3; i++) {
            anchor_coord(0, i) = _scaffold(m1[anchor[0].num], i);
            anchor_coord(1, i) = _scaffold(m1[anchor[2].num], i);
            anchor_coord(2, i) = _scaffold(m1[anchor[3].num], i);
            anchor_coord(3, i) = _scaffold(m1[anchor[1].num], i);
        }
        int len = anchor[2].num - anchor[0].num + 1;
        auto helix_coord = make_helix(anchor_coord, len);
        std::cout << "helix:\n" << helix_coord << std::endl;
        MatrixXf temp_mat = mat::hstack(helix_coord.row(0), helix_coord.row(len - 1));
        temp_mat = mat::hstack(temp_mat, helix_coord.row(len));
        temp_mat = mat::hstack(temp_mat, helix_coord.row(2 * len - 1));
        sp(helix_coord, temp_mat, anchor_coord);
        std::cout << "superposed helix:\n" << helix_coord << std::endl;
        helix_coords[m1[anchor[0].num] / 2] = helix_coord.topRows(helix_coord.rows() / 2);
        helix_coords[m1[anchor[3].num] / 2] = helix_coord.bottomRows(helix_coord.rows() / 2);
    }
//    std::cout << "\nCalculate helix coordinate: " << std::endl;
//    std::vector<MatrixXf> helix_coords;
//    SupPos sp;
//    for (int i = 0; i < _scaffold.rows(); i += 2) {
//        auto helix_coord = make_helix_strand(_scaffold.row(i), _scaffold.row(i + 1), m2[i + 1] - m2[i]);
//        MatrixXf m, n;
//        m = mat::hstack(helix_coord.row(0), helix_coord.row(helix_coord.rows() - 1));
//        n = mat::hstack(_scaffold.row(i), _scaffold.row(i + 1));
//        sp(helix_coord, m, n);
//        helix_coords.push_back(helix_coord);
//        std::cout << "helix coordinate: \n" << helix_coord << std::endl;
//    }

    /// Construct fragments
    std::vector<int> anchor_list;
    for (auto &&vec: _helix_anchors) {
        for (auto &&res: vec) {
            anchor_list.push_back(res.num - 1);
        }
    }
    std::sort(anchor_list.begin(), anchor_list.end(), [](int a, int b){return a < b;});
    auto fragments = get_fragments(anchor_list);

    /// Calculate coordinates of fragments
    MatrixXf coords;
    int frag_num = 0;
    for (auto &&frag: fragments) {
        MatrixXf a, b;
        if (std::get<0>(frag).size() != 0) {
            a.resize(2, 3);
            for (int i = 0; i < 3; i++) {
                a(0, i) = _scaffold(m1[std::get<0>(frag)[0] + 1], i);
                a(1, i) = _scaffold(m1[std::get<0>(frag)[1] + 1], i);
            }
        }
        if (std::get<2>(frag).size() != 0) {
            b.resize(2, 3);
            for (int i = 0; i < 3; i++) {
                b(0, i) = _scaffold(m1[std::get<2>(frag)[0] + 1], i);
                b(1, i) = _scaffold(m1[std::get<2>(frag)[1] + 1], i);
            }
        }
        auto frag_coords = get_frag_coords(frag, a, b);
        int temp_len = a.rows() + b.rows();
        MatrixXf x(temp_len, 3), y(temp_len, 3);
        if (std::get<0>(frag).size() != 0) {
            for (int i = 0; i < 3; i++) {
                x(0, i) = frag_coords(0, i);
                x(1, i) = frag_coords(1, i);
                y(0, i) = a(0, i);
                y(1, i) = a(1, i);
            }
        }
        if (std::get<2>(frag).size() != 0) {
            for (int i = 0; i < 3; i++) {
                x(x.rows() - 2, i) = frag_coords(frag_coords.rows() - 2, i);
                x(x.rows() - 1, i) = frag_coords(frag_coords.rows() - 1, i);
                y(y.rows() - 2, i) = b(0, i);
                y(y.rows() - 1, i) = b(1, i);
            }
        }
        sp(frag_coords, x, y);
        std::cout << "after superpose: " << std::endl;
        std::cout << frag_coords << std::endl;
//        if (std::get<2>(frag).size() != 0) {
//            coords = mat::hstack(coords, frag_coords.topRows(frag_coords.rows() - 2));
//        } else {
//            coords = mat::hstack(coords, frag_coords);
//        }
        if (std::get<0>(frag)[0] == 0) {
            coords = mat::hstack(helix_coords[frag_num], frag_coords.block(2, 0, frag_coords.rows() - 4, 3));
            frag_num++;
            coords = mat::hstack(coords, helix_coords[frag_num]);
        } else if (frag_num == 0) {
            coords = mat::hstack(coords, frag_coords.topRows(frag_coords.rows() - 2));
            coords = mat::hstack(coords, helix_coords[frag_num]);
        } else if (frag_num == fragments.size() - 1) {
            coords = mat::hstack(coords, frag_coords.bottomRows(frag_coords.rows() - 2));
        } else {
            coords = mat::hstack(coords, frag_coords.block(2, 0, frag_coords.rows() - 4, 3));
            coords = mat::hstack(coords, helix_coords[frag_num]);
        }

        frag_num++;
    }


    std::cout << "\nfinal coords: \n" << coords << std::endl;

}

MatrixXf LM2::make_helix(const MatrixXf &anchor, int len) {
    int num_atoms = len * 2;
    MatrixXf bound(num_atoms, num_atoms);
    for (int i = 0; i < num_atoms; i++) {
        for (int j = i; j < num_atoms; j++) {
            if (i == j) {
                bound(i, j) = 0;
            } else if (i < len && j < len || i >= len && j >= len) {
                int d = j - i;
                bound(i, j) = bound(j, i) = sqrt(2*9.7*9.7*(1-cos(0.562*d-2*3.14159))+(2.84*d)*(2.84*d));
            } else if (i < len && j >= len && j < num_atoms - 1 - i) {
                int d = num_atoms - 1 - i - j;
                bound(i, j) = bound(j, i) = sqrt(2*9.7*9.7*(1-cos(0.562*d-1.5*3.14159))+(2.84*d-4)*(2.84*d-4));
            } else if (i < len && j >= len && j > num_atoms - 1 - i) {
                int d = j - (num_atoms - 1 - i);
                bound(i, j) = bound(j, i) = sqrt(2*9.7*9.7*(1-cos(0.562*d-0.5*3.14159))+(2.84*d+4)*(2.84*d+4));
            } else if (j == num_atoms - 1 - i) {
                bound(i, j) = bound(j, i) = 15.1;
            }
        }
    }
    std::cout << "helix bound:\n" << bound << std::endl;
    DG dg(bound);
    auto coord = dg();
    std::cout << "helix energy: " << dg.E << std::endl;
    return coord;
}

MatrixXf LM2::make_helix_strand(MatrixXf head, MatrixXf tail, int len) {
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

MatrixXf LM2::get_frag_coords(const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> &frag, MatrixXf a, MatrixXf b) {
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
            } else if (m2[j] - m2[i] == 1) {
                bound(i, j) = bound(j, i) = 6.1;
            } else {
                bound(i, j) = 6.1 * (m2[j] - m2[i]);
                bound(j, i) = 6.1;
            }
        }
    }
    if (a.rows() == 0) {
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
    std::cout << "\nfragment bound: \n" << bound << std::endl;
    DG dg(bound);
    auto coord = dg();
//    std::cout << coord << "\n" << "energy: " << dg.E << "\n" << std::endl;
    std::cout << "DG...\nenergy: " << dg.E << "\n" << coord << std::endl;
    std::cout << "distances: " << std::endl;
    for (int i = 0; i < coord.rows() - 1; i++) {
        std::cout << (coord.row(i) - coord.row(i + 1)).norm() << ' ';
    }
    std::cout << std::endl;
    return coord;
}

std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> LM2::get_fragments(const std::vector<int> &anchor_list) {
    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> fragments;
    if (anchor_list[0] != 0) {
        std::vector<int> a, b, c;
        for (int j = 0; j < anchor_list[0]; j++) {
            b.push_back(j);
        }
        c.push_back(anchor_list[0]);
        c.push_back(anchor_list[1]);
        fragments.push_back(std::make_tuple(a, b, c));
    }
    for (int i = 0; i < anchor_list.size() / 2 - 1; i++) {
        std::vector<int> a, b, c;
        a.push_back(anchor_list[i * 2]);
        a.push_back(anchor_list[i * 2 + 1]);
        for (int j = anchor_list[i * 2 + 1] + 1; j < anchor_list[i * 2 + 2]; j++) {
            b.push_back(j);
        }
        c.push_back(anchor_list[i * 2 + 2]);
        c.push_back(anchor_list[i * 2 + 3]);
        fragments.push_back(std::make_tuple(a, b, c));
    }
    if (anchor_list.back() != _seq.size() - 1) {
        std::vector<int> a, b, c;
        a.push_back(anchor_list[anchor_list.size() - 2]);
        a.push_back(anchor_list.back());
        for (int j = anchor_list.back() + 1; j < _seq.size(); j++) {
            b.push_back(j);
        }
        fragments.push_back(std::make_tuple(a, b, c));
    }

    std::cout << "\nfragments: " << std::endl;
    for (auto &&frag: fragments) {
        for (auto &&i: std::get<0>(frag)) {
            std::cout << i << '-';
        }
        for (auto &&i: std::get<1>(frag)) {
            std::cout << i << '-';
        }
        for (auto &&i: std::get<2>(frag)) {
            std::cout << i << '-';
        }
        std::cout << ' ';
        std::cout << std::endl;
    }
    return fragments;
}


void LM2::set_helix_coords(loop *src) {
    for (bp *b = src->s.head; b != NULL; b = b->next) {
        if (std::count_if(_helix_anchors.begin(), _helix_anchors.end(), [&](const std::vector<res> &vec){return std::count_if(vec.begin(), vec.end(), [&](const res &r2){return b->res1.type == r2.type && b->res1.name == r2.name && b->res1.num == r2.num;});})) {
            std::cout << b->res1.num << ' ' << b->res1.type << ' ' << b->res1.name << std::endl;
        }
        if (std::count_if(_helix_anchors.begin(), _helix_anchors.end(), [&](const std::vector<res> &vec){return std::count_if(vec.begin(), vec.end(), [&](const res &r2){return b->res2.type == r2.type && b->res2.name == r2.name && b->res2.num == r2.num;});})) {
            std::cout << b->res2.num << ' ' << b->res2.type << ' ' << b->res2.name << std::endl;
        }
    }
    for (res *r = src->head; r != NULL; r = r->next) {
        if (std::count_if(_helix_anchors.begin(), _helix_anchors.end(), [&](const std::vector<res> &vec){return std::count_if(vec.begin(), vec.end(), [&](const res &r2){return r->type == r2.type && r->name == r2.name && r->num == r2.num;});})) {
            std::cout << r->num << ' ' << r->type << ' ' << r->name << std::endl;
        }
    }
}

void LM2::set_helix_anchors(loop *src) {
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

R5P LM2::get_helix(const helix &h) {
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

R5P LM2::create_helix(const string &seq) {
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
        return RNA(file_name);
    } else {
        string file_name = _lib + "/basepair/" + seq.substr(0, 3) + 
                           seq.substr(seq.size() - 3, 3) + ".pdb";
        ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = _lib + "/basepair/XXXXXX.pdb";
        }
        ifile.close();
        return Connect()(RNA(file_name), create_helix(seq.substr(1, seq.size() - 2)), 2, 3);
    }
}

MatrixXf LM2::get_helix_par(const R5P &r5p) {
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

void LM2::set_constraints() {
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
        map<int, int> temp_map;
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

MatrixXf LM2::apply_constraints(const MatrixXf &scaffold) {
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
        map<int, int> temp_map;
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

std::vector<std::vector<std::vector<int>>> LM2::read_constraints_file() {
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

MatrixXf LM2::read_constraints_pos(std::string type, int layer_nums) {
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

MatrixXf LM2::read_constraints_par(std::string type, int layer_nums) {
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

} /// namespace nuc3d
} /// namespace jian


