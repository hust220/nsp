#ifndef JIAN_PDB_MANIPULATE_H
#define JIAN_PDB_MANIPULATE_H

#include "Pdb.h"
#include "../geom/geometry.h"
#include "../geom/rotate.h"

namespace jian {
    
namespace pdb {
    
inline Vector3f normal_vector(const Residue &res) {
    if (std::set<std::string>{"A", "G", "DA", "DG"}.count(res.name)) {
        auto atom_list = res[std::vector<std::string>{"C2", "C4", "C6"}];
        return geometry::normal_vector(atom_list[0], atom_list[1], atom_list[2]);
    } else if (std::set<std::string>{"U", "C", "DT", "DC"}.count(res.name)) {
        auto atom_list = res[std::vector<std::string>{"C2", "C6", "C4"}];
        return geometry::normal_vector(atom_list[0], atom_list[1], atom_list[2]);
    }
}

inline Vector3f pivot(const Residue &res) {
    if (std::set<std::string>{"A", "G", "DA", "DG"}.count(res.name)) {
        return res["N9"].pos<Vector3f>();
    } else if (std::set<std::string>{"U", "C", "DT", "DC"}.count(res.name)) {
        return res["N1"].pos<Vector3f>();
    }
}
 
} // namespace pdb

} // namespace jian



#endif

