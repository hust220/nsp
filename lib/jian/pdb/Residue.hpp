#pragma once

#include <string>
#include <deque>
#include "../utils/traits.hpp"
#include "Atom.hpp"

BEGIN_JN

	class Residue : public std::deque<Atom> {
	public:
		int num = -1;
		S name = "X";
		S m_cg = "aa";

		JN_DEFAULT_CONSTRUCTORS(Residue);

		static auto get_sort_keys();

		void sort();

		S format_name(const S &s);

		Atom &operator [](int n);

		const Atom &operator [](int n) const;

		Atom &operator [](const S &s);

		const Atom &operator [](const S &s) const;

		Residue &operator +=(const Residue &r);

		bool has_atom(const Str &atom_name) const {
			return std::find_if(this->begin(), this->end(), [&atom_name](const Atom &atom) {
				return atom.name == atom_name;
			}) != this->end();
		}

		template<typename _Second, typename... _Rest>
		bool has_atom(const Str &first, _Second &&second, _Rest &&...rest) const {
			 return has_atom(first) && has_atom(second, rest...);
		}

		template<typename _Atoms>
		void set_atoms(const _Atoms &atoms) {
			this->clear();
			for (const Atom &atom : atoms) {
				this->push_back(atom);
			}
		}

		JN_DEF_ATOMS;
	};

#define JN_DEF_RESIDUES \
	refs<Residue> residues() { return refs<Residue>().append(*this); }\
	refs<const Residue> residues() const { return refs<const Residue>().append(*this); }


	bool res_is_type(const Residue &res, S type);

	template<typename T>
	uniform_const_t<Atom, T> &atom(T &&res, const S &s) {
		for (auto &&atom : res) if (atom.name == s) return atom;
		throw "jian::pdb::Atom::operator [](const S &) error!";
	}

	template<typename T, template<typename...> class L = std::deque>
	auto atoms(const Residue &r, T &&ls) {
		L<Atom> v;
		for (auto &&atom : r) if (std::count(ls.begin(), ls.end(), atom.name)) v.push_back(atom);
		return v;
	}

	template<typename T>
	bool exists_atom(T &&res, const S &s) {
		for (auto &&atom : res) if (atom.name == s) return true;
		return false;
	}

END_JN

