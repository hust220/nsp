#pragma once

#include "pdb_molecule.hpp"
#include "jian.hpp"

namespace jian {

#define JN_IS_ATOM(T)     std::is_same<std::decay_t<T>, Atom    >::value
#define JN_IS_RESIDUE(T)  std::is_same<std::decay_t<T>, Residue >::value
#define JN_IS_CHAIN(T)    std::is_same<std::decay_t<T>, Chain   >::value
#define JN_IS_MODEL(T)    std::is_same<std::decay_t<T>, Model   >::value
#define JN_IS_MOLECULE(T) std::is_same<std::decay_t<T>, Molecule>::value

#define JN_GE_ATOM(T)     JN_IS_ATOM(T)     || JN_GE_RESIDUE(T)
#define JN_GE_RESIDUE(T)  JN_IS_RESIDUE(T)  || JN_GE_CHAIN(T)
#define JN_GE_CHAIN(T)    JN_IS_CHAIN(T)    || JN_GE_MODEL(T)
#define JN_GE_MODEL(T)    JN_IS_MODEL(T)    || JN_GE_MOLECULE(T)
#define JN_GE_MOLECULE(T) JN_IS_MOLECULE(T)

inline void chain_append(Chain &chain, const Chain &c) {
	for (auto && r : c) {
		chain.push_back(r);
	}
}

inline void chain_append(Chain &chain, Chain &&c) {
	for (auto && r : c) {
		chain.push_back(std::move(r));
	}
}

inline void chain_append(Chain &chain, const Model &model) {
	for (auto && c : model) {
		for (auto && r : c) {
			chain.push_back(r);
		}
	}
}

inline void chain_append(Chain &chain, Model &&model) {
	for (auto && c : model) {
		for (auto && r : c) {
			chain.push_back(std::move(r));
		}
	}
}

// seq
template<typename T>
Str seq(T &&t) {
	Str s;
	each<Residue>(t, [&s](auto &&res) {
		s += res.name;
		return true;
	});
	return s;
}

template<typename T>
int num_residues(T &&t) { return count<Residue>(t); }

template<typename T>
int num_atoms(T &&t) { return count<Atom>(t); }

template<typename T, JN_ENABLE(JN_GE_RESIDUE(T))>
int is_empty(T && t) {
	return num_atoms(t) == 0;
}

// sub
template<typename K>
Molecule sub(const Molecule &mol, K &&k) {
	Molecule m;
	for (auto && model : mol) {
		Model && mod = sub(model, std::forward<K>(k));
		if (!mod.empty()) {
			m.push_back(std::move(mod));
		}
	}
	return m;
}

template<typename K>
Model sub(const Model & model, K && t) {
	int res_num = 0;
	Model m;

	for (auto && chain : model) {
		Chain temp_chain;
		temp_chain.name = chain.name;
		for (auto &&res : chain) {
			if (std::find(std::begin(t), std::end(t), res_num) != std::end(t)) {
				temp_chain.push_back(res);
			}
			res_num++;
		}
		m.push_back(temp_chain);
	}

	return m;
}

template<typename K>
Chain sub(const Chain & chain, K && t) {
	int res_num = 0;
	Chain temp_chain;

	temp_chain.name = chain.name;
	for (auto &&res : chain) {
		if (std::find(std::begin(t), std::end(t), res_num) != std::end(t)) {
			temp_chain.push_back(res);
		}
		res_num++;
	}

	return temp_chain;
}

template<template<typename...> class L = std::deque>
auto residues(const Model &model) {
	L<Residue> v;
	for (auto &&chain : model) for (auto &&res : chain) v.push_back(res);
	return v;
}

template<typename T, template<typename...> class L = std::deque>
auto residues(const Model &model, T &&ls) {
	L<Residue> v;
	int i = 0; for (auto &&chain : model) for (auto &&res : chain) {
		if (std::count(std::begin(ls), std::end(ls), i)) v.push_back(res);
		i++;
	}
	return v;
}

template<typename T>
uniform_const_t<Residue, T> &residue(T &&model, int n) {
	int i = 0; for (auto &&chain : model) for (auto &&res : chain) {
		if (i == n) return res;
		i++;
	}
	throw "JIAN::MODEL::residue(int) error! Residue index out of range.";
}

template<typename T>
auto residues_to_model(T &&ls) {
	Model model; model.resize(1);
	for (auto &&res : ls) model[0].push_back(res);
	return model;
}

}


