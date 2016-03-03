#ifndef CHAIN_H_INCLUDED
#define CHAIN_H_INCLUDED

#include "Residue.h"

namespace jian {

class Chain : public std::deque<Residue> {
public:
    std::string name = "X";
    std::string type = "unknown";

    Chain() {}

    Chain(MolFile &pdb_file) {
        if (!pdb_file.eof()) {
            name = pdb_file.chain_name();
            int model_num = pdb_file.model_num();
            while (!pdb_file.eof() && model_num == pdb_file.model_num() && name == pdb_file.chain_name()) {
                this->push_back(Residue(pdb_file));
            }
        }
    }

};

} // namespace jian

#endif 

