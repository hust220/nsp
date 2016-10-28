#pragma once

#include "Molecule.hpp"
#include "MolParser.hpp"
#include "../utils/file.hpp"

namespace jian {

void chain_read_model(Chain &chain, std::string f, std::string type = "");
Chain read_model_to_chain(std::string f, std::string type = "");
void append_chain_to_file(const Chain &chain, const std::string &file_name, int n);

template<typename T>
void mol_write(const T & t, std::string file_name) {
    std::ofstream(file_name.c_str()) << t;
}

template<typename T>
void mol_read(T &t, std::string file_name, std::string mol_type = "") {
    MolParser *parser = MolParser::make(file::type(file_name), file_name, mol_type);
    (*parser) >> t;
    delete parser;
}

template<typename T>
T mol_read_to(std::string f, std::string type = "") {
    T t;
    mol_read(t, f, type);
    return t;
}

MolParser &operator >>(MolParser &input, Atom &atom);
MolParser &operator >>(MolParser &input, Residue &residue);
MolParser &operator >>(MolParser &input, Chain &chain);
MolParser &operator >>(MolParser &input, Model &model);
MolParser &operator >>(MolParser &input, Molecule &mol);

std::ostream &operator <<(std::ostream &output, const Atom &atom);
std::ostream &operator <<(std::ostream &output, const Residue &residue);
std::ostream &operator <<(std::ostream &output, const Chain &chain);
std::ostream &operator <<(std::ostream &output, const Model &model);
std::ostream &operator <<(std::ostream &output, const Molecule &mol);

template<typename T, typename U, typename = std::enable_if_t<!(std::is_same<T, Residue>::value), int>>
T coarse_grained(const T & t, const U & ls) {
    T t_;
    t_.name = t.name;
    for (auto && r : t) {
        t_.push_back(coarse_grained(r, ls));
    }
    return t_;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_same<T, Residue>::value, int>>
T coarse_grained(const T &residue, U && ls) {
    Residue r;
    r.name = residue.name;
    r.num = residue.num;
    for (auto && atom : residue) {
        if (std::find(ls.begin(), ls.end(), atom.name) != ls.end()) {
            r.push_back(atom);
        }
    }
    return r;
}

} // namespace jian

