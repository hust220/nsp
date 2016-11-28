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
		std::string m_cg = "aa";

		JN_DEFAULT_CONSTRUCTORS(Residue);

		static auto get_sort_keys();

		void sort();

		std::string format_name(const std::string &s);

		Atom &operator [](int n);

		const Atom &operator [](int n) const;

		Atom &operator [](const std::string &s);

		const Atom &operator [](const std::string &s) const;

		Residue &operator +=(const Residue &r);

		JN_DEF_ATOMS;
	};

#define JN_DEF_RESIDUES \
	refs<Residue> residues() { return refs<Residue>().append(*this); }\
	refs<const Residue> residues() const { return refs<const Residue>().append(*this); }


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

