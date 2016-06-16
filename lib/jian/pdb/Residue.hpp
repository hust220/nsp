#pragma once

#include <string>
#include <deque>
#include "../utils/traits.hpp"
#include "Atom.hpp"

namespace jian {

class Residue : public std::deque<Atom> {
public:
    int num = -1;
    std::string name = "X";
    std::string type = "unknown";

    Residue() {}
    static auto get_sort_keys();
    void sort();
    std::string format_name(const std::string &s);
    Atom &operator [](int n);
    const Atom &operator [](int n) const;
    Atom &operator [](const std::string &s);
    const Atom &operator [](const std::string &s) const;

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

Residue residue_from_file(const std::string &f);

void residue_to_file(const Residue &residue, const std::string &file_name);

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

