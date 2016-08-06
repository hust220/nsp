#pragma once

#include <string>
#include <vector>
#include <map>

namespace jian {
namespace pdb {

using names_t = std::vector<std::string>;
using map_names_t = std::map<std::string, names_t>;

struct names {
    static names_t res;
    static names_t bases;
    static map_names_t atoms_res;
    static names_t atoms_phos;
    static map_names_t atoms_sugar;
    static map_names_t atoms_base;
};

}
}

