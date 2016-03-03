#ifndef JIAN_PDB_CONVERT_H
#define JIAN_PDB_CONVERT_H

#include "Pdb.h"
#include "../geom.h"
#include "Format.h"

namespace jian {

class Convert {
public:
    using Mat = MatrixXd;
    std::string _lib;
    std::map<std::string, std::vector<std::string>> _atom_list;

    Convert() {
        char *lib = getenv("NSP");
        if (lib == NULL) die("Please set environment variable 'NSP'!");
        _lib += lib;
        _atom_list["A"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", 
                           "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"};
        _atom_list["U"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", 
                           "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"};
        _atom_list["G"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", 
                           "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"};
        _atom_list["C"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", 
                           "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"};
        _atom_list["DA"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", 
                            "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"};
        _atom_list["DT"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", 
                            "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"};
        _atom_list["DG"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", 
                            "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"};
        _atom_list["DC"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", 
                            "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"};
    }

    template<class T> Atom make_atom(std::string name, const T &coord) {
        Atom atom;
        atom.name = name;
        for (int i = 0; i < 3; i++) atom[i] = coord[i];
        return atom;
    }

    template<class T> Residue make_residue(std::string name, const T &coord) {
        Residue res;
        res.name = name;
        int index = (_atom_list[name].size() == coord.rows() ? 0 : 3);
        for (int i = 0; i < coord.rows(); i++) {
            res.push_back(
                make_atom(_atom_list[name][index], coord.row(i)));
            index++;
        }
        return res;
    }

    void operator ()(Residue &res, std::string name) {
        if (res.name == name) return;

        /// Sort atoms of res
        Format format;
        format.sort(res);

        /// Get information
        Mat res_coord; 
        std::vector<int> res_nail; 
        std::vector<int> res_direct;
        std::tie(res_coord, res_nail, res_direct) = get_res_coord(res);

        Mat base_coord; 
        std::vector<int> base_nail; 
        std::vector<int> base_direct;
        std::tie(base_coord, base_nail, base_direct) = read_base(name);
     
        /// Superpose axis
        geom::suppos(base_coord, base_nail, res_coord, res_nail);

        /// Calculate base direction
        RowVector3d direct1 = geom::normal_vector(
            base_coord.row(base_direct[0]),
            base_coord.row(base_direct[1]),
            base_coord.row(base_direct[2])
        );
        RowVector3d direct2 = geom::normal_vector(
            res_coord.row(res_direct[0]),
            res_coord.row(res_direct[1]),
            res_coord.row(res_direct[2])
        );


        /// Calculate rotate angle
        auto rotate_angle = geom::dihedral(
            base_coord.row(base_nail[0]) + direct1,
            base_coord.row(base_nail[0]),
            base_coord.row(base_nail[1]),
            base_coord.row(base_nail[1]) + direct2
        );

        /// Rotate base
        VectorXd beg = base_coord.row(base_nail[0]);
        VectorXd end = base_coord.row(base_nail[1]);
        geom::rotate_angle_xyz(base_coord, beg, end, rotate_angle);

        /// Set residue coordiantes size
        std::set<std::string> dna_name_list{"DA", "DT", "DG", "DC"};
        std::set<std::string> rna_name_list{"A", "U", "G", "C"};

        int old_phos_size, old_sugar_size, old_base_size;
        std::tie(old_phos_size, old_sugar_size, old_base_size) = get_res_size(res);

        int new_phos_size, new_sugar_size, new_base_size;
        new_phos_size = old_phos_size;
        new_sugar_size = (dna_name_list.count(name) ? 8 : 9);
        BOOST_ASSERT(new_sugar_size == old_sugar_size && "jian::Convert error! Please check the residue to be converted.");
        new_base_size = std::map<std::string, int>{{"A", 10}, {"U", 8}, {"G", 11}, {"C", 8}, {"DA", 10}, {"DT", 9}, {"DG", 11}, {"DC", 8}}[name];

        Mat coord(new_phos_size + new_sugar_size + new_base_size, 3);

        /// Copy phosphate and sugar coordinates
        if (rna_name_list.count(res.name) && dna_name_list.count(name)) {
            coord.topRows(new_phos_size + new_sugar_size - 1) = res_coord.topRows(old_phos_size + old_sugar_size - 2);
        } else if (dna_name_list.count(res.name) && rna_name_list.count(name)) {
            coord.topRows(new_phos_size + new_sugar_size - 2) = res_coord.topRows(old_phos_size + old_sugar_size - 1);
            coord.row(new_phos_size + new_sugar_size - 2) = get_o2_coord(res);
        } else {
            coord.topRows(new_phos_size + new_sugar_size) = res_coord.topRows(old_phos_size + old_sugar_size);
        }

        /// Copy base coordinates
        coord.bottomRows(base_coord.rows()) = base_coord;

        /// Make residue
        res = make_residue(name, coord);
    }

    std::tuple<Mat, std::vector<int>, std::vector<int>> get_res_coord(const Residue &res) {
        Mat res_coord(res.size(), 3);
        std::vector<int> res_nail(2);
        std::vector<int> res_direct(3);
        for (int i = 0; i < res.size(); i++) {
            const Atom &atom = res[i];
            for (int j = 0; j < 3; j++) res_coord(i, j) = atom[j];
            if (atom.name == "C1*") {
                res_nail[0] = i;
            } else if (std::set<std::string>{"A", "G", "DA", "DG"}.count(res.name) && atom.name == "N9") {
                res_nail[1] = i;
            } else if (std::set<std::string>{"U", "C", "DT", "DC"}.count(res.name) && atom.name == "N1") {
                res_nail[1] = i;
            } else if (atom.name == "C2") {
                res_direct[0] = i;
            } else if (atom.name == "C4") {
                res_direct[1] = i;
            } else if (atom.name == "C6") {
                res_direct[2] = i;
            }
        }
        if (std::set<std::string>{"U", "C", "DT", "DC"}.count(res.name)) {
            std::swap(res_direct[1], res_direct[2]);
        }
        return std::make_tuple(res_coord, res_nail, res_direct);
    }

    void set_res_coord(Residue &res, const Mat &res_coord) {
        for (int i = 0; i < res.size(); i++) for (int j = 0; j < 3; j++) res[i][j] = res_coord(i, j);
    }

    std::tuple<Mat, std::vector<int>, std::vector<int>> read_base(std::string name) {
        std::string file_name = _lib;
        if (std::set<std::string>{"A", "U", "G", "C"}.count(name)) {
            file_name += "/RNA/base/" + name + ".pdb";
        } else if (std::set<std::string>{"DA", "DT", "DG", "DC"}.count(name)) {
            file_name += "/DNA/base/" + name + ".pdb";
        }
        return get_res_coord(Model(file_name)[0][0]);
    }

    std::tuple<int, int, int> get_res_size(const Residue &res) {
        std::set<std::string> phos_atoms{"P", "O1P", "O2P"};
        int phos_size = 0, sugar_size = 0, base_size = 0;
        for (auto &&atom: res) {
            if (phos_atoms.count(atom.name)) {
                phos_size++;
            } else if (atom.name.size() == 3) {
                sugar_size++;
            } else {
                base_size++;
            }
        }
        return make_tuple(phos_size, sugar_size, base_size);
    }

    RowVector3d get_o2_coord(const Residue &res) {
        RowVector3d c1, c2, c3, o2;
        for (auto &&atom: res) {
            if (atom.name == "C1*") { for (int i = 0; i < 3; i++) c1[i] = atom[i];
            } else if (atom.name == "C2*") { for (int i = 0; i < 3; i++) c2[i] = atom[i];
            } else if (atom.name == "C3*") { for (int i = 0; i < 3; i++) c3[i] = atom[i];
            }
        }
        o2 = (c3 - c2).cross(c1 - c2) + (c3 + c1 - 2 * c2) + c2;
        return o2;
    }

};

}

#endif


