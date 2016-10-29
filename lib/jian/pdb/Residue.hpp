#pragma once

#include <string>
#include <deque>
#include "../utils/traits.hpp"
#include "Atom.hpp"

namespace jian {

class Residue : public std::deque<Atom> {
public:
	enum cg_code {
		CG_AA,
		CG_1P,
		CG_6P,
		CG_PSB
	};

    int num;
    std::string name;

	Residue();
    static auto get_sort_keys();
    void sort();
    std::string format_name(const std::string &s);
    Atom &operator [](int n);
    const Atom &operator [](int n) const;
    Atom &operator [](const std::string &s);
    const Atom &operator [](const std::string &s) const;

	template<typename CG_T>
	bool is_cg() const {
		return m_cg == CG_T::m_cg;
	}

	template<typename CG_T>
	Residue & cg() {
		*this = CG_T::res(*this);
		m_cg = CG_T::m_cg;
		return *this;
	}

	cg_code m_cg;
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

