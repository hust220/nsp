#include "Format.h"

namespace jian {

namespace pdb {

Residue Format::operator ()(const Residue &res) {
    Residue new_res;
    new_res.name = res.name;
    std::vector<std::string> v{"P", "O1P", "O2P"};
    if (std::all_of(v.begin(), v.end(), [&](std::string s){return std::any_of(res.atoms.begin(), res.atoms.end(), [&](const Atom &atom){return atom.name == s;});})) {
        new_res.atoms.push_back(res["P"]);
        new_res.atoms.push_back(res["O1P"]);
        new_res.atoms.push_back(res["O2P"]);
    }
    for (auto &&name_pair: _atom_rank[res.name]) {
        if (std::count(v.begin(), v.end(), name_pair.first)) continue;
        if (std::none_of(res.atoms.begin(), res.atoms.end(), [&](const Atom &atom){return atom.name == name_pair.first;})) {
            return Residue();
        } else {
            new_res.atoms.push_back(res[name_pair.first]);
        }
    }
    sort(new_res);
    return new_res;
//    if (std::all_of(_atom_rank[res.name].begin(), _atom_rank[res.name].end(), [&](const std::pair<std::string, int> name_pair){return std::set<std::string>{"P", "O1P", "O2P"}.count(name_pair.first) || std::any_of(res.atoms.begin(), res.atoms.end(), [&](const Atom &atom){return atom.name == name_pair.first;});})) {
//        new_res = res;
//        sort(new_res);
//    }
//    return new_res;
}

Chain Format::operator ()(const Chain &chain) {
    Chain new_chain;
    new_chain.name = chain.name;
    for (auto &&res: chain.residues){
        auto r = (*this)(res);
        if (! r.atoms.empty()) {
            new_chain.residues.push_back(r);
        }
    }   
    return new_chain;
}

Model Format::operator ()(const Model &model) {
    Model new_model;
    new_model.name = model.name;
    for (auto &&chain: model.chains) {
        auto c = (*this)(chain);
        if (! c.residues.empty()) {
            new_model.chains.push_back(c);
        }
    }
    return new_model;
}

Pdb Format::operator ()(const Pdb &pdb) {
    Pdb new_pdb;
    new_pdb.name = pdb.name;
    for (auto &&model: pdb.models) {
        auto m = (*this)(model);
        if (! m.chains.empty()) {
            new_pdb.models.push_back(m);
        }
    }
    return new_pdb;
}

void Format::sort(Residue &res) {
    std::sort(res.atoms.begin(), res.atoms.end(), [&](const Atom &atom1, const Atom &atom2) {
        return _atom_rank[res.name][atom1.name] < _atom_rank[res.name][atom2.name];
    });
}

void Format::sort(Chain &chain) {
    for (auto &&residue: chain.residues) {
        (*this)(residue);
    }
}

void Format::sort(Model &mol) {
    for (auto &&chain: mol.chains) {
        (*this)(chain);
    }
}

void Format::sort(Pdb &pdb) {
    for (auto &&mol: pdb.models) {
        (*this)(mol);
    }
}

} /// namespace pdb

} /// namespace jian


