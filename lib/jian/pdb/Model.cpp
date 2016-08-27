#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <map>
#include "Model.hpp"
#include "molstream.hpp"
#include "../utils/file.hpp"

namespace jian {

std::string seq(const Model &model) {
    std::string seq;
    for (auto &&chain : model) for (auto &&res : chain) { seq += res.name; }
    return seq;
}

int num_residues(const Model &model) {
    int i = 0; for (auto &&chain : model) for (auto &&res : chain) i++;
    return i;
}

int num_atoms(const Model &model) {
    int i = 0; for (auto &&chain : model) for (auto &&res : chain) for (auto &&atom : res) i++;
    return i;
}

bool is_empty(const Model &model) {
    return num_residues(model) == 0;
}

//Model RNA(const Model &model) {
//    Model rna;
//    thread_local static std::set<std::string> names {"A", "U", "G", "C"};
//    for (auto &&chain: model) {
//        Chain temp_chain; temp_chain.name = chain.name;
//        for (auto &&residue: chain) {
//            auto res = residue; res.name = res.name.substr(0, 1);
//            if (names.count(res.name)) temp_chain.push_back(std::move(res));
//        }
//        if (!temp_chain.empty()) rna.push_back(temp_chain);
//    }
//    rna.name = model.name; rna.type = "RNA";
//    return rna;
//}
//
//Model RNA(const std::string &s) {
//    return read_model(s, "RNA");
//}
//
//Model DNA(const Model &model) {
//    Model dna;
//    thread_local static std::set<std::string> names {"DA", "DT", "DG", "DC"};
//    for (auto &&chain: model) {
//        Chain temp_chain; temp_chain.name = chain.name;
//        for (auto &&residue: chain) {
//            auto res = residue; res.name = res.name.substr(0, 2);
//            if (names.count(res.name)) temp_chain.push_back(std::move(res));
//        }
//        if (!temp_chain.empty()) dna.push_back(temp_chain);
//    }
//    dna.name = model.name; dna.type = "DNA";
//    return dna;
//}
//
//Model DNA(const std::string &s) {
//    return DNA(Model(s));
//}
//
//Model R5P(const Model &model) {
//    thread_local static std::map<std::string, std::set<std::string>> names {
//        {"A", {"C5*", "O3*", "C1*", "N6", "C2"}},
//        {"U", {"C5*", "O3*", "C1*", "O2", "O4"}},
//        {"G", {"C5*", "O3*", "C1*", "O6", "N2"}},
//        {"C", {"C5*", "O3*", "C1*", "O2", "N4"}}
//    };
//    Model m; m.type = "R5P";
//    for (auto &&chain : model) {
//        Chain new_chain; new_chain.name = chain.name; new_chain.type = m.type;
//        for (auto &&res : chain) {
//            Residue new_residue; new_residue.name = res.name;
//            for (auto &&atom : res) if (names[res.name].count(atom.name)) new_residue.push_back(atom);
//            if (!new_residue.empty()) new_chain.push_back(new_residue);
//        }
//        if (!new_chain.empty()) m.push_back(new_chain);
//    }
//    return m;
//}
//
//Model R5P(const std::string &s) {
//    return R5P(RNA(s));
//}

} // namespace jian

