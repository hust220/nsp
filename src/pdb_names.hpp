#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace jian {
namespace pdb {

using names_t = std::vector<std::string>;
using map_names_t = std::map<std::string, names_t>;

class Names {
public:
    using map_Names = std::map<std::string, std::shared_ptr<Names>>;

    static const Names & instance(S mol_type);

    names_t res;
    map_names_t alias;
    map_names_t atoms_base;
    map_names_t atoms_res;
    names_t atoms_phos;
    names_t atoms_sugar;
    names_t atoms_bb;

private:
    Names();
    static void init_map(map_Names &map);
    void print_names(const Names & names);
};

inline int res_type(const S &rname) {
    for (auto && mol_type : { "RNA", "DNA", "protein" }) {
        const Names &names = Names::instance(mol_type);
        auto it = std::find_if(names.res.begin(), names.res.end(), [&rname, &names](const S &s) {
            const names_t &v = names.alias.at(s);
            return s == rname || std::find(v.begin(), v.end(), rname) != v.end();
        });
        if (it != names.res.end()) {
            return std::distance(names.res.begin(), it);
        }
    }
    throw std::string("unknown residue: ") + rname;
}


inline Str res_name(const Str &rname, Str type = "std") {
    for (auto && mol_type : { "RNA", "DNA", "protein" }) {
        const Names &names = Names::instance(mol_type);
        auto it = std::find_if(names.res.begin(), names.res.end(), [&rname, &names](const S &s) {
            const names_t &v = names.alias.at(s);
            return s == rname || std::find(v.begin(), v.end(), rname) != v.end();
        });
        if (it != names.res.end()) {
            if (it->size() == 1 || type == "std") return *it;
            else {
                const names_t &v = names.alias.at(*it);
                it = std::find_if(v.begin(), v.end(), [](auto && s){return s.size() == 1;});
                return *it;
            }
        }
    }
    return "X";
}

inline int res_mol_type(const S &rname) {
    int i = 0;
    for (auto && mol_type : { "RNA", "DNA", "protein" }) {
        const Names &names = Names::instance(mol_type);
        auto it = std::find_if(names.res.begin(), names.res.end(), [&rname, &names](const S &s) {
            const names_t &v = names.alias.at(s);
            return s == rname || std::find(v.begin(), v.end(), rname) != v.end();
        });
        if (it != names.res.end()) {
            return i;
        }
        i++;
    }
    throw std::string("unknown residue: ") + rname;
}

template<typename T>
inline T res_mol_type(const S &rname) {
    for (auto && mol_type : { "RNA", "DNA", "protein" }) {
        //std::cout << mol_type << std::endl;
        const Names &names = Names::instance(mol_type);
        auto it = std::find_if(names.res.begin(), names.res.end(), [&rname, &names](const S &s) {
            const names_t &v = names.alias.at(s);
            return s == rname || std::find(v.begin(), v.end(), rname) != v.end();
        });
        if (it != names.res.end()) {
            return mol_type;
        }
    }
    throw std::string("unknown residue: ") + rname;
}

inline const names_t &res_included_atoms(const S &rname) {
    return Names::instance(res_mol_type<std::string>(rname)).atoms_res.at(res_name(rname, "abbr"));
}

}

}

