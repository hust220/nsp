#ifndef JIAN_NUC3D_BUILDSTRAND_H
#define JIAN_NUC3D_BUILDSTRAND_H

#include "../pdb/Residue.h"
#include "../util/std.h"
#include "../util/mat.h"
#include "../dg/DG.h"
#include "../geom/geometry.h"

namespace jian {
namespace nuc3d {

class BuildStrand {
public:
    Log log;
    std::set<std::string> _coarse_atoms {"C4*"};
    DG dg;
    SupPos sp;
    std::string _lib = env("NSP");
    std::string _frag_par_path;
    std::vector<std::string> _frag_3_names;
    std::vector<std::string> _frag_4_names;
    std::vector<std::string> _frag_5_names;
    std::vector<std::array<double,  3>> _frag_3_dists;
    std::vector<std::array<double,  6>> _frag_4_dists;
    std::vector<std::array<double, 10>> _frag_5_dists;


    BuildStrand() {
        read_pars();
        dg.log.set_display(false);
    }

    std::vector<Residue> operator ()(int n, const MatrixXf &a, const MatrixXf &b) {
        return build_strand(n, a, b);
    }

    std::vector<Residue> build_strand(int n, const MatrixXf &a, const MatrixXf &b) {
        if ((a.rows() != _coarse_atoms.size() * 2 && a.rows() != 0) || (b.rows() != _coarse_atoms.size() * 2 && b.rows() != 0)) {
            throw "jian::nuc3d::BuildStrand::build_strand(int, const MatrixXf &, const MatrixXf &) error!";
        }
        log("build strand:\n");
        auto bound = make_bound(n, a, b);
        log("matrix a:\n", a, "\n", "matrix b:\n", b, "\n", "bound:\n", bound, "\n");
        auto scaffold = dg(bound);
        log("scaffold:\n", scaffold, "\n");
        superpose_scaffold(scaffold, a, b);
        log("scaffold after superposed:\n", scaffold, "\n", scaffold.rows(), '\n');
        auto residues = all_atom(scaffold);
        auto new_residues = slice(residues, a.rows(), a.rows() + n);
        if (a.rows() != 0) {
            for (int i = 0; i < 3; i++) log(new_residues.front()["C4*"][i], ' ');
            log('\n');
        }
        if (b.rows() != 0) {
            for (int i = 0; i < 3; i++) log(new_residues.back()["C4*"][i], ' ');
            log('\n');
        }
        return new_residues;
        //return slice(residues, a.rows(), a.rows() + n);
        //return slice(all_atom(scaffold), a.rows(), a.rows() + n);
    }

    MatrixXf make_bound(int n, const MatrixXf &a, const MatrixXf &b) {
        int size = n * _coarse_atoms.size() + a.rows() + b.rows();
        MatrixXf bound(size, size);
        bound = MatrixXf::Zero(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                if (i == j) {
                    bound(i, j) = 0;    
                } else {
                    bound(i, j) = 6.1 * (j - i);
                    bound(j, i) = 6.1;
                }
            }
        }
        auto mat = mat::hstack(a, b);
        std::vector<int> vec;
        if (a.rows() != 0 && b.rows() != 0) vec = std::vector<int> {0, 1, size - 2, size - 1};
        else if (a.rows() == 0) vec = std::vector<int> {size - 2, size - 1};
        else if (b.rows() == 0) vec = std::vector<int> {0, 1};
        else return bound;
        for (int i = 0; i < mat.rows(); i++) {
            for (int j = i + 1; j < mat.rows(); j++) {
                bound(vec[i], vec[j]) = bound(vec[j], vec[i]) = geometry::distance(mat.row(i), mat.row(j));
            }
        }
        return bound;
    }

    void superpose_scaffold(MatrixXf &scaffold, const MatrixXf &a, const MatrixXf &b) {
        if (a.rows() != 0 || b.rows() != 0) {
            int temp_len = a.rows() + b.rows();
            MatrixXf x(temp_len, 3), y(temp_len, 3);
            if (a.rows() != 0) {
                for (int i = 0; i < 3; i++) {
                    x(0, i) = scaffold(0, i);
                    x(1, i) = scaffold(1, i);
                    y(0, i) = a(0, i);
                    y(1, i) = a(1, i);
                }
            }
            if (b.rows() != 0) {
                for (int i = 0; i < 3; i++) {
                    x(x.rows() - 2, i) = scaffold(scaffold.rows() - 2, i);
                    x(x.rows() - 1, i) = scaffold(scaffold.rows() - 1, i);
                    y(y.rows() - 2, i) = b(0, i);
                    y(y.rows() - 1, i) = b(1, i);
                }
            }
            sp(scaffold, x, y);
        }
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
            log(_frag_par_path, "/", _frag_3_names[rank.top().first], ".pdb\n");
            return get_residues_from_file(_frag_par_path + "/" + _frag_3_names[rank.top().first] + ".pdb");
        } else if (size == 10) {
            for (int i = 0; i < _frag_5_dists.size(); i++) 
                rank.push(std::make_pair(i, geometry::norm(_frag_5_dists[i], std::forward<List>(list), 10)));
            log(_frag_par_path, "/", _frag_5_names[rank.top().first], ".pdb\n");
            return get_residues_from_file(_frag_par_path + "/" + _frag_5_names[rank.top().first] + ".pdb");
        } else {
            throw "JIAN::NUC3D::LM2::find_best_frag_model error!";
        }
    }

    template<typename Coord> std::vector<double> coord_to_dists(Coord &&coord) {
        int len = coord.rows();
        std::vector<double> dists(len * (len - 1) / 2);
        for (int i = 0, index = 0; i < len; i++) for (int j = i + 1; j < len; j++) {
            dists[index] = geometry::distance(coord.row(i), coord.row(j));
            index++;
        }
        return dists;
    }

    template<typename Residues, typename Coord> void superpose_residues(Residues &&residues, Coord &&coord) {
        int len = coord.rows();
        MatrixXf x(len, 3), y(len, 3);
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
            x(i, j) = residues[i]["C4*"][j];
            y(i, j) = coord(i, j);
        }
        sp(x, y);
        auto c1 = -sp.c1;
        for (auto &&res: residues) for (auto &&atom: res.atoms) {
            geom::move(atom, c1);    
            geom::rotate(atom, sp.rot);    
            geom::move(atom, sp.c2);    
        }
    }

    std::vector<Residue> all_atom(const MatrixXf &scaffold) {
        int len = scaffold.rows();
        std::vector<Residue> residues;
        residues.reserve(len);
        if (len >= 5) {
            for (int i = 0; i < len - 4; i++) {
                auto coord = scaffold.block(i, 0, 5, 3);
                auto frag_residues = find_best_frag_model(coord_to_dists(coord));
                superpose_residues(frag_residues, coord);
                if (i == 0) append(residues, frag_residues[0], frag_residues[1]);
                residues.push_back(frag_residues[2]);
                if (i == len - 5) append(residues, frag_residues[3], frag_residues[4]);
            }
            return residues;
        } else {
            auto frag_residues = find_best_frag_model(coord_to_dists(scaffold));
            superpose_residues(frag_residues, scaffold);
            for (auto &&res: frag_residues) residues.push_back(std::move(res));
            return residues;
        }
    }

    void read_pars() {
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

};

} // namespace nuc3d
} // namespace jian

#endif

