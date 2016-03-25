#pragma once

#include "Atom.hpp"

namespace jian {

class Residue : public std::deque<Atom> {
public:
    int num = -1;
    std::string name = "X";
    std::string type = "unknown";

    Residue() {}

    Residue(MolFile &pdb_file) {
        if (!pdb_file.eof()) {
            name = format_name(pdb_file.res_name());
            num = pdb_file.res_num();
            int model_num = pdb_file.model_num();
            std::string chain_name = pdb_file.chain_name();
            while (!pdb_file.eof() and model_num == pdb_file.model_num() and chain_name == pdb_file.chain_name() 
                   and num == pdb_file.res_num() and name == format_name(pdb_file.res_name())) {
                this->push_back(Atom(pdb_file));
            }
        }
    }

    static auto get_sort_keys() {
        static std::map<std::string, std::vector<std::string>> keys {
            {"A",{"P","O1P","O2P",
                  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
                  "N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"}},
            {"U",{"P","O1P","O2P",
                  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
                  "N1","C2","O2","N3","C4","O4","C5","C6"}},
            {"G",{"P","O1P","O2P",
                  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
                  "N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"}},
            {"C",{"P","O1P","O2P",
                  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
                  "N1","C2","O2","N3","C4","N4","C5","C6"}}
        };
        std::map<std::string, std::map<std::string, int>> sort_keys;
        int index = 0;
        for (auto && res_name : {"A", "U", "G", "C"}) {
            for (int i = 0; i < keys[res_name].size(); i++) {
                sort_keys[res_name][keys[res_name][i]] = index;
                index++;
            }
        }
        return sort_keys;
    }

    void sort() {
        static auto sort_keys = get_sort_keys();
        auto & keys = sort_keys;
        std::sort(this->begin(), this->end(), [&](auto &&a1, auto &&a2){
            return keys[name][a1.name] < keys[name][a2.name];
        });
    }

                  
    std::string format_name(const std::string &s) {
        std::smatch result;
        if (std::regex_match(s, result, std::regex("^(\\w+)\\d+$"))) {
            // tLeap would append a '5' after the name of first residue
            // and '3' after the name of the last residue
            return result[1];
        } else return s;
    }              
                   
    Atom &operator [](int n) {
        return std::deque<Atom>::operator [](n);
    }              
                   
    const Atom &operator [](int n) const {
        return std::deque<Atom>::operator [](n);
    }

    Atom &operator [](const std::string &s) {
        for (auto &&atom : *this) if (atom.name == s) {return atom;}
        throw "jian::Residue::operator[] error! Not found atom!";
    }

    const Atom &operator [](const std::string &s) const {
        for (auto &&atom : *this) if (atom.name == s) {return atom;}
        throw "jian::Residue::operator[] error! Not found atom!";
    }

    template<typename T>
    Residue coarse_grained(T &&names) const {
        Residue r;
        r.name = name;
        for (auto &&atom : *this) {
            if (std::find(names.begin(), names.end(), atom.name) != names.end()) {
                r.push_back(atom);
            }
        }
        return r;
    }
};

template<typename T>
uniform_const_t<Atom, T> &atom(T &&res, const std::string &s) {
    for (auto &&atom : res) if (atom.name == s) return atom;
    throw "jian::pdb::Atom::operator [](const std::string &) error!";
}

template<typename T, template<typename...> class L = std::deque>
auto atoms(const Residue &r, T &&ls) {
    L<Atom> v;
    for (auto &&atom : r) if (std::count(ls.begin(), ls.end(), atom.name)) v.push_back(atom);
    return v;
}

template<typename T>
bool exists_atom(T &&res, const std::string &s) {
    for (auto &&atom : res) if (atom.name == s) return true;
    return false;
}

} // namespace jian

