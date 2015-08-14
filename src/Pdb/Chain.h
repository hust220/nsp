#ifndef CHAIN_H_INCLUDED
#define CHAIN_H_INCLUDED

#include "Residue.h"

namespace jian {

class Chain
{
public:
    Chain();
    Chain(MolFile &mol_file);
//    Chain(PdbFile &pdb_file);
//    Chain(Cif &cif);
    Chain(Chain *chain) : name(chain->name), residues(chain->residues) {}
    Chain(const Chain &chain) : name(chain.name), residues(chain.residues) {}
    Chain &operator =(const Chain &chain) {
        name = chain.name;
        residues = chain.residues;
        return *this;
    }
    Chain(vector<string> &lines,  string type = "");
    int atom_nums();

    void push(Residue *);
    void push(Residue &);
    Residue &operator [](int);
    const Residue &operator [](int) const;

    string name = "A";
    vector<Residue> residues;

    vector<Residue>::iterator begin();
    vector<Residue>::iterator end();
};

} /// namespace jian

#endif // CHAIN_H_INCLUDED
