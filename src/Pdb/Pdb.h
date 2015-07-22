#ifndef Pdb_H
#define Pdb_H

#include "Model.h"

namespace jian {

class Pdb
{
public:
    Pdb() {}
    Pdb(string pdbfile) {
        readPDB(pdbfile);
    }

    vector<Model>::iterator begin() {
        return models.begin();
    }
    vector<Model>::iterator end() {
        return models.end();
    }

    void push(const Model &);
    void readPDB(string);
    Model &operator [](int);
    const Model &operator [](int) const;

    /* attributes */
    string name = "none";
    vector<Model> models;

};

} /// namespace jian

#endif // Pdb_H

