#include <algorithm>
#include <regex>
#include "../utils/file.hpp"
#include "../utils/Env.hpp"
#include "names.hpp"

namespace jian {
namespace pdb {


class Names {
public:
    names_t bases = {"A", "U", "T", "G", "C"};
    names_t res;
    names_t atoms_phos;
    map_names_t atoms_res, atoms_sugar, atoms_base;

    Names() {
        std::string file_name = Env::lib() + "/RNA/pars/pdb/names";
        std::vector<std::string> v;
        EACH_LINE(file_name.c_str(),
            if (std::regex_search(L, std::regex("Phos"))) {
                std::getline(ifile, L);
                jian::tokenize(L, atoms_phos, " ");
            } else if (std::regex_search(L, std::regex("Sugar"))) {
                if (std::regex_search(L, std::regex("RNA"))) {
                    std::getline(ifile, L);
                    jian::tokenize(L, atoms_sugar["RNA"], " ");
                } else if (std::regex_search(L, std::regex("DNA"))) {
                    std::getline(ifile, L);
                    jian::tokenize(L, atoms_sugar["DNA"], " ");
                }
            } else if (std::regex_search(L, std::regex("Base\\s+A"))) {
                std::getline(ifile, L);
                jian::tokenize(L, atoms_base["A"], " ");
            } else if (std::regex_search(L, std::regex("Base\\s+T"))) {
                std::getline(ifile, L);
                jian::tokenize(L, atoms_base["T"], " ");
            } else if (std::regex_search(L, std::regex("Base\\s+U"))) {
                std::getline(ifile, L);
                jian::tokenize(L, atoms_base["U"], " ");
            } else if (std::regex_search(L, std::regex("Base\\s+G"))) {
                std::getline(ifile, L);
                jian::tokenize(L, atoms_base["G"], " ");
            } else if (std::regex_search(L, std::regex("Base\\s+C"))) {
                std::getline(ifile, L);
                jian::tokenize(L, atoms_base["C"], " ");
            }
        );
        for (auto && name : bases) {
            res.push_back(name);
            std::copy(atoms_phos.begin(), atoms_phos.end(), std::back_inserter(atoms_res[name]));
            std::copy(atoms_sugar["RNA"].begin(), atoms_sugar["RNA"].end(), std::back_inserter(atoms_res[name]));
            std::copy(atoms_base[name].begin(), atoms_base[name].end(), std::back_inserter(atoms_res[name]));
            std::string s = std::string("D") + name;
            res.push_back(s);
            std::copy(atoms_phos.begin(), atoms_phos.end(), std::back_inserter(atoms_res[s]));
            std::copy(atoms_sugar["DNA"].begin(), atoms_sugar["DNA"].end(), std::back_inserter(atoms_res[s]));
            std::copy(atoms_base[name].begin(), atoms_base[name].end(), std::back_inserter(atoms_res[s]));
        }
    }
};

thread_local static Names impl;

names_t names::bases = impl.bases;
names_t names::res = impl.res;
map_names_t names::atoms_res = impl.atoms_res;
names_t names::atoms_phos = impl.atoms_phos;
map_names_t names::atoms_sugar = impl.atoms_sugar;
map_names_t names::atoms_base = impl.atoms_base;

} // namespace pdb
} // namespace jian

