#pragma once

#include "../pdb.hpp"
#include <vector>
#include <map>
#include <string>

namespace jian {

class Convert {
public:
    using Mat = Eigen::MatrixXd;
    std::string _lib;
    std::map<std::string, std::vector<std::string>> _atom_list;

    Convert();
    void operator ()(Residue &res, std::string name);

    void set_res_coord(Residue &res, const Mat &res_coord);
    auto get_res_coord(const Residue &res);
    auto read_base(std::string name);
    auto get_res_size(const Residue &res);
    auto get_o2_coord(const Residue &res);

    template<class T> 
    Atom make_atom(std::string name, const T &coord) {
        Atom atom;
        atom.name = name;
        for (int i = 0; i < 3; i++) atom[i] = coord[i];
        return atom;
    }

    template<class T> 
    Residue make_residue(std::string name, const T &coord) {
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

};

} // namespace jian

