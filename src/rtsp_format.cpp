#include <vector>
#include <algorithm>
#include "rtsp_format.hpp"

BEGIN_JN

Format::Format() {
    _atom_rank["A"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"O2*", 10}, {"C1*", 11},
                      {"N9", 12}, {"C8", 13}, {"N7", 14}, {"C5", 15},
                      {"C6", 16}, {"N6", 17}, {"N1", 18}, {"C2", 19},
                      {"N3", 20}, {"C4", 21}};
    _atom_rank["U"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"O2*", 10}, {"C1*", 11},
                      {"N1", 12}, {"C2", 13}, {"O2", 14}, {"N3", 15},
                      {"C4", 16}, {"O4", 17}, {"C5", 18}, {"C6", 19}};
    _atom_rank["G"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"O2*", 10}, {"C1*", 11},
                      {"N9", 12}, {"C8", 13}, {"N7", 14}, {"C5", 15},
                      {"C6", 16}, {"O6", 17}, {"N1", 18}, {"C2", 19},
                      {"N2", 20}, {"N3", 21}, {"C4", 22}};
    _atom_rank["C"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"O2*", 10}, {"C1*", 11},
                      {"N1", 12}, {"C2", 13}, {"O2", 14}, {"N3", 15},
                      {"C4", 16}, {"N4", 17}, {"C5", 18}, {"C6", 19}};
    _atom_rank["DA"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"C1*", 10},
                      {"N9", 11}, {"C8", 12}, {"N7", 13}, {"C5", 14},
                      {"C6", 15}, {"N6", 16}, {"N1", 17}, {"C2", 18},
                      {"N3", 19}, {"C4", 20}};
    _atom_rank["DT"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"C1*", 10},
                      {"N1", 11}, {"C2", 12}, {"O2", 13}, {"N3", 14},
                      {"C4", 15}, {"O4", 16}, {"C5", 17}, {"C7", 18},
                      {"C6", 19}};
    _atom_rank["DG"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"C1*", 10},
                      {"N9", 11}, {"C8", 12}, {"N7", 13}, {"C5", 14},
                      {"C6", 15}, {"O6", 16}, {"N1", 17}, {"C2", 18},
                      {"N2", 19}, {"N3", 20}, {"C4", 21}};
    _atom_rank["DC"] = {{"P", 0}, {"O1P", 1}, {"O2P", 2}, {"O5*", 3},
                      {"C5*", 4}, {"C4*", 5}, {"O4*", 6}, {"C3*", 7},
                      {"O3*", 8}, {"C2*", 9}, {"C1*", 10},
                      {"N1", 11}, {"C2", 12}, {"O2", 13}, {"N3", 14},
                      {"C4", 15}, {"N4", 16}, {"C5", 17}, {"C6", 18}};
}

Residue Format::operator ()(const Residue &res) {
    Str name_std = pdb::res_name(res.name, "std");
    Residue new_res;
    new_res.name = name_std;
    std::vector<std::string> v{"P", "O1P", "O2P"};
    if (std::all_of(v.begin(), v.end(), [&](S s){
		return std::any_of(res.begin(), res.end(), [&](const Atom &atom){
			return atom.name == s;
		});
	})) {
        new_res.push_back(atom(res, "P"));
        new_res.push_back(atom(res, "O1P"));
        new_res.push_back(atom(res, "O2P"));
    }
    for (auto &&name_pair: _atom_rank[name_std]) {
        if (std::count(v.begin(), v.end(), name_pair.first)) continue;
        if (std::none_of(res.begin(), res.end(), [&](const Atom &atom){
			return atom.name == name_pair.first;
		})) {
            return Residue();
        } else {
            new_res.push_back(atom(res, name_pair.first));
        }
    }
    sort(new_res);
    return new_res;
}

Chain Format::operator ()(const Chain &chain) {
    Chain new_chain;
    new_chain.name = chain.name;
    for (auto &&res: chain){
        auto r = (*this)(res);
        if (! r.empty()) {
            new_chain.push_back(r);
        }
    }   
    return new_chain;
}

Model Format::operator ()(const Model &model) {
    Model new_model;
    new_model.name = model.name;
    for (auto &&chain: model) {
        auto c = (*this)(chain);
        if (! c.empty()) {
            new_model.push_back(c);
        }
    }
    return new_model;
}

Molecule Format::operator ()(const Molecule &pdb) {
    Molecule new_pdb;
    new_pdb.name = pdb.name;
    for (auto &&model : pdb) {
        auto m = (*this)(model);
        if (! m.empty()) {
            new_pdb.push_back(m);
        }
    }
    return new_pdb;
}

void Format::sort(Residue &res) {
    std::sort(res.begin(), res.end(), [&](const Atom &atom1, const Atom &atom2) {
        return _atom_rank[res.name][atom1.name] < _atom_rank[res.name][atom2.name];
    });
}

void Format::sort(Chain &chain) {
    for (auto &&residue: chain) {
        this->sort(residue);
    }
}

void Format::sort(Model &mol) {
    for (auto &&chain: mol) {
        this->sort(chain);
    }
}

void Format::sort(Molecule &pdb) {
    for (auto &&mol: pdb) {
        this->sort(mol);
    }
}

END_JN

