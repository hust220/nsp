#pragma once

#include <string>
#include <vector>
#include <map>

namespace jian {
namespace pdb {

using names_t = std::vector<std::string>;
using map_names_t = std::map<std::string, names_t>;

class Names {
public:
	Names();
	static Names & instance();

    names_t res;
    names_t bases;
    map_names_t atoms_res;
    names_t atoms_phos;
    map_names_t atoms_sugar;
    map_names_t atoms_base;
};

}
}

