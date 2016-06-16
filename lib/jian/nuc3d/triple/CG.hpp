#pragma once

#include <string>
#include <vector>
#include <map>

namespace jian {
namespace triple {

struct CG {
    using atom_t = std::string;
    using res_t = std::string;
    using atoms_t = std::vector<atom_t>;
    using res_atoms_map_t = std::map<res_t, atoms_t>;

    static res_atoms_map_t res_atoms_map;
};

} // namespace triple
} // namespace jian

