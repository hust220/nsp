#ifndef Pdb_H
#define Pdb_H

#include "Model.h"

namespace jian {

class Pdb {
public:
    Pdb();
    Pdb(MolFile &mol_file);
//    Pdb(PdbFile &pdb_file);
//    Pdb(Cif &cif);
    Pdb(string file_name);

    vector<Model>::iterator begin() {
        return models.begin();
    }
    vector<Model>::iterator end() {
        return models.end();
    }

    void push(const Model &);
    void read(string);
    void read_pdb(string);
    Model &operator [](int);
    const Model &operator [](int) const;

    /* attributes */
    string name = "none";
    vector<Model> models;

};

} /// namespace jian

#endif // Pdb_H

