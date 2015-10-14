#ifndef DNA_H
#define DNA_H

#include "Model.h"

namespace jian {


class DNA :public Model {
public:
    DNA();
    DNA(const Model &model);
    DNA(DNA *dna) {
        name = dna->name;
        chains = dna->chains;
    }
    DNA(const DNA &dna) {
        name = dna.name;
        chains = dna.chains;
    }
    DNA &operator =(const DNA &dna) {
        name = dna.name;
        chains = dna.chains;
        return *this;
    }
    friend ostream &operator <<(ostream &, const DNA &);
    void print();
    void write(string);
    DNA(char *);
    DNA(string);
    DNA *copy();
    void read_pdb(string);
    void read_cif(string);
    void push(Chain *);
    void push(Chain &);
    Chain &operator [](int);
    void updateChains(string);

    void setResNum();

    int getLen();
    int totalAtoms();
    double getDist(int, int);
    string getSeq();
    string getChain();
    
    string name;
};

} /// namespace jian

#endif // DNA_H

