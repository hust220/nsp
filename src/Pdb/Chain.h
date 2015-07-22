#ifndef CHAIN_H_INCLUDED
#define CHAIN_H_INCLUDED

#include "Residue.h"

namespace jian {

class Chain
{
public:
    Chain() {}
    Chain(Chain *chain) : name(chain->name), rnaName(chain->rnaName), residues(chain->residues) {}
    Chain(const Chain &chain) : name(chain.name), rnaName(chain.rnaName), residues(chain.residues) {}
    Chain &operator =(const Chain &chain) {
        name = chain.name;
        rnaName = chain.rnaName;
        residues = chain.residues;
        return *this;
    }
    Chain(vector<string> &, string = "");
    int atom_nums();

    void push(Residue *);
    void push(Residue &);
    Residue &operator [](int);
    const Residue &operator [](int) const;

    string name = "A";
    string rnaName;
    vector<Residue> residues;

    vector<Residue>::iterator begin();
    vector<Residue>::iterator end();
};

} /// namespace jian

#endif // CHAIN_H_INCLUDED
