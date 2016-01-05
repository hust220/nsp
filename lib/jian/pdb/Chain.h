#ifndef CHAIN_H_INCLUDED
#define CHAIN_H_INCLUDED

#include "Residue.h"

namespace jian {

template<typename ResType>
class BasicChain {
public:
    using ResidueType = ResType;

    std::string name = "A";
    std::vector<ResType> residues;

    BasicChain() {}
    BasicChain(BasicChain<ResType> *chain) : name(chain->name), residues(chain->residues) {}
    BasicChain(const BasicChain<ResType> &chain) : name(chain.name), residues(chain.residues) {}

    BasicChain<ResType> &operator =(const BasicChain<ResType> &chain) {
        name = chain.name;
        residues = chain.residues;
        return *this;
    }

    BasicChain(MolFile &pdb_file) {
        if (!pdb_file.eof()) {
            name = pdb_file.chain_name();
            int model_num = pdb_file.model_num();
            while (!pdb_file.eof() && model_num == pdb_file.model_num() && name == pdb_file.chain_name()) {
                residues.push_back(ResType(pdb_file));
            }
        }
    }

    BasicChain(std::vector<std::string> &lines, std::string type = "RNA") {
        /// Set name
        name = "";
        name += lines[0][21];

        /// set residues
        vector<string> strings;
        for (int i = 0; i < (int) lines.size(); i++) {
            if (!strings.empty() && lines[i].substr(22, 5) != strings.back().substr(22, 5)) {
                ResType res(strings);
                if (boost::to_upper_copy(type) == "RNA") {
                    if (set<string>{"A", "U", "G", "C"}.count(res.name)) {
                        residues.push_back(ResType(strings));
                    }
                } else if (boost::to_upper_copy(type) == "DNA") {
                    if (set<string>{"DA", "DT", "DG", "DC"}.count(res.name)) {
                        residues.push_back(ResType(strings));
                    }
                } else {
                    residues.push_back(ResType(strings));
                }
                strings.clear();
            }
            strings.push_back(lines[i]);
        }
        ResType res(strings);
        if (boost::to_upper_copy(type) == "RNA") {
            if (set<string>{"A", "U", "G", "C"}.count(res.name)) {
                residues.push_back(ResType(strings));
            }
        } else if (boost::to_upper_copy(type) == "DNA") {
            if (set<string>{"DA", "DT", "DG", "DC"}.count(res.name)) {
                residues.push_back(ResType(strings));
            }
        } else {
            residues.push_back(ResType(strings));
        }
        strings.clear();
    }

    int atom_nums() {
        int num = 0;
        for (auto &&residue: residues) {
            for (auto &&atom: residue) {
                num++;
            }
        }
        return num;
    }

    void push(ResType *residue) {
        residues.push_back(*residue);
    }

    void push(ResType &residue) {
        residues.push_back(residue);
    }

    ResType &operator [](int n) {
        return residues[n];
    }

    const ResType &operator [](int n) const {
        return residues[n];
    }

    typename std::vector<ResType>::iterator begin() {
        return residues.begin();
    }

    typename std::vector<ResType>::iterator end() {
        return residues.end();
    }
};

typedef BasicChain<Residue> Chain;

} // namespace jian

#endif 

