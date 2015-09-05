#ifndef CONVERT_H
#define CONVERT_H

#include "Pdb.h"
#include "../Geom.h"
#include "Format.h"

namespace jian {

class Convert {
public:
    Convert() {
        char *lib = getenv("NSP");
        if (lib == NULL) die("Please set environment variable 'NSP'!");
        _lib += lib;
        _atom_list["A"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"};
        _atom_list["U"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"};
        _atom_list["G"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"};
        _atom_list["C"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"};
        _atom_list["DA"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"};
        _atom_list["DT"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"};
        _atom_list["DG"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"};
        _atom_list["DC"] = {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "C1*", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"};
    }

    void operator ()(Residue &, std::string);
    std::tuple<MatrixXf, std::vector<int>, std::vector<int>> get_res_coord(const Residue &);
    void set_res_coord(Residue &, const MatrixXf &);
    std::tuple<MatrixXf, std::vector<int>, std::vector<int>> read_base(std::string);
    std::tuple<int, int, int> get_res_size(const Residue &);
    RowVector3f get_o2_coord(const Residue &);

    template<class T>
    Atom make_atom(std::string name, const T &coord) {
        Atom atom;
        atom.name = name;
        atom.x = coord[0];
        atom.y = coord[1];
        atom.z = coord[2];
        return atom;
    }

    template<class T>
    Residue make_residue(std::string name, const T &coord) {
        Residue res;
        res.name = name;
        int index = 
            (_atom_list[name].size() == coord.rows() ? 0 : 3);
        for (int i = 0; i < coord.rows(); i++) {
            res.atoms.push_back(
                make_atom(_atom_list[name][index], coord.row(i)));
            index++;
        }
        return res;
    }

private:
    std::string _lib;
    std::map<std::string, std::vector<std::string>> _atom_list;
};

}

#endif


