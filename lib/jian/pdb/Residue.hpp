#pragma once

#include <string>
#include <deque>
#include "../utils/traits.hpp"
#include "Atom.hpp"

namespace jian {

class Residue : public std::deque<Atom> {
public:
    int num;
    std::string name;
	std::string m_cg;

	Residue();
    static auto get_sort_keys();
    void sort();
    std::string format_name(const std::string &s);
    Atom &operator [](int n);
    const Atom &operator [](int n) const;
    Atom &operator [](const std::string &s);
    const Atom &operator [](const std::string &s) const;

};

bool res_is_type(const Residue &res, std::string type);

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

