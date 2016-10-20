#pragma once

#include <map>
#include <string>
#include "../pdb.hpp"

namespace jian {

class Format {
private:
    std::map<std::string, std::map<std::string, int>> _atom_rank;

public:
    Format();
    Residue operator ()(const Residue &res);
    Chain operator ()(const Chain &chain);
    Model operator ()(const Model &model);
    Molecule operator ()(const Molecule &pdb);
    void sort(Residue &res);
    void sort(Chain &chain);
    void sort(Model &mol);
    void sort(Molecule &pdb);
};

} // namespace jian

