#ifndef RNA_H
#define RNA_H

#include "Model.h"

namespace jian {

class RNA :public Model {
public:
    RNA() : Model() {}
    RNA(const Model &model): Model(model) {}
    RNA(RNA *rna) {
        len = rna->len;
        chains = rna->chains;
    }
    RNA(const RNA &rna) {
        len = rna.len;
        chains = rna.chains;
    }
    RNA &operator =(const RNA &rna) {
        len = rna.len;
        chains = rna.chains;
        return *this;
    }
    RNA(char *);
    RNA(string);
    RNA *copy();
    void read(string);
    void read_pdb(string);
    void read_cif(string);
    void push(Chain *);
    void push(Chain &);
    Chain &operator [](int);
    void updateChains(string);

    void setLen();
    void setResNum();

    int getLen();
    int totalAtoms();
    double getDist(int, int);
    string getSeq();
    string getChain();
    
    /* IO function */
    void print();
    void printAsDNA();
    void write(string);
    
    /* assemble function */
    void move(double, double, double);
    void rotate(Matr_ *m);
    void format();
    void addP();
    void mutate(string);
    void rotateByX(double);
    void rotateByZ(double);

    /* attributes */
//    string name;
    int len;
};

class RNAs {
public:
    RNAs();
    ~RNAs();

    int getLen();
    void resize(int);
    RNA *at(int);
    RNA *operator [](int);
    void push(RNA *);
    
    vector<RNA *> RNAList;
};

} /// namespace jian

#endif // RNA_H
