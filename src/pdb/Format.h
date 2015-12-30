#ifndef JIAN_FORMAT_H
#define JIAN_FORMAT_H

#include "Pdb.h"

namespace jian {

namespace pdb {

class Format {
private:
    std::map<std::string, std::map<std::string, int>> _atom_rank;

public:
    Format() {
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

    Residue operator ()(const Residue &res) {
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
    }

    Chain operator ()(const Chain &chain) {
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

    Model operator ()(const Model &model) {
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

    Pdb operator ()(const Pdb &pdb) {
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

    void sort(Residue &res) {
        std::sort(res.atoms.begin(), res.atoms.end(), [&](const Atom &atom1, const Atom &atom2) {
            return _atom_rank[res.name][atom1.name] < _atom_rank[res.name][atom2.name];
        });
    }

    void sort(Chain &chain) {
        for (auto &&residue: chain.residues) {
            (*this)(residue);
        }
    }

    void sort(Model &mol) {
        for (auto &&chain: mol.chains) {
            (*this)(chain);
        }
    }

    void sort(Pdb &pdb) {
        for (auto &&mol: pdb.models) {
            (*this)(mol);
        }
    }

};

} /// namespace pdb

} /// namespace jian

#endif

