#ifndef JIAN_PDB_MANIPULATE_H
#define JIAN_PDB_MANIPULATE_H

#include "Pdb.hpp"
#include "../geom.hpp"

namespace jian {
    
inline Vector3d normal_vector(const Residue &res) {
    if (std::set<std::string>{"A", "G", "DA", "DG"}.count(res.name)) {
        auto atom_list = atoms(res, std::vector<std::string>{"C2", "C4", "C6"});
        return geom::normal_vector(atom_list[0], atom_list[1], atom_list[2]);
    } else if (std::set<std::string>{"U", "C", "DT", "DC"}.count(res.name)) {
        auto atom_list = atoms(res, std::vector<std::string>{"C2", "C6", "C4"});
        return geom::normal_vector(atom_list[0], atom_list[1], atom_list[2]);
    }
}

inline Vector3d pivot(const Residue &res) {
    if (std::set<std::string>{"A", "G", "DA", "DG"}.count(res.name)) {
        return pos<Vector3d>(atom(res, "N9"));
    } else if (std::set<std::string>{"U", "C", "DT", "DC"}.count(res.name)) {
        return pos<Vector3d>(atom(res, "N1"));
    }
}
 
} // namespace jian



#endif

