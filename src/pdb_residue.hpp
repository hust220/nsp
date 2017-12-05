#pragma once

#include <string>
#include <deque>
#include "jian.hpp"
#include "pdb_atom.hpp"

BEGIN_JN

class Residue : public std::deque<Atom> {
public:
	int num = -1;
	S name = "X";
	S m_cg = "aa";

	JN_DEFAULT_CONSTRUCTORS(Residue);

	Atom &operator [](int n);

	const Atom &operator [](int n) const;

	Atom &operator [](const S &s);

	const Atom &operator [](const S &s) const;

	Residue &operator +=(const Residue &r);

	JN_DEF_ATOMS;
};

#define JN_DEF_RESIDUES \
	refs<Residue> residues() { return refs<Residue>().append(*this); }\
	refs<const Residue> residues() const { return refs<const Residue>().append(*this); }

Str format_res_name(const Str &s);

bool is_mol_type(const Residue &res, Str type);

void sort(Residue &res);

template<typename _Ls>
void set_atoms(Residue &res, const _Ls &atoms) {
	res.clear();
	for (const Atom &atom : atoms) {
		res.push_back(atom);
	}
}

inline Bool has_atom(const Residue &res, const Str &atom_name) {
	return std::find_if(res.begin(), res.end(), [&atom_name](const Atom &atom) {
		return atom.name == atom_name;
	}) != res.end();
}

template<typename... _Rest>
Bool has_atom(const Residue &res, const Str &first, const Str &second, _Rest &&...rest) {
	return has_atom(res, first) && has_atom(res, second, rest...);
}

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

END_JN

