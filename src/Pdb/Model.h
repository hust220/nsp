#ifndef Model_H
#define Model_H

#include "Chain.h"

namespace jian {

class Model
{
public:
    Model() {}
    Model(Model *model) : name(model->name), chains(model->chains) {}
    Model(const Model &model) : name(model.name), chains(model.chains) {}
    Model &operator =(const Model &model) {
        name = model.name;
        chains = model.chains;
    }
    Model(vector<string>);
    Model(std::string pdbfile);

    vector<Chain>::iterator begin();
    vector<Chain>::iterator end();

    int res_nums() const {
        int res_num = 0;
        for (auto &chain: chains) {
            for (auto &residue: chain.residues) {
                res_num++;
            }
        }
        return res_num;
    }

    int atom_nums() const {
        int atom_num = 0;
        for (auto &chain: chains) {
            for (auto &residue: chain.residues) {
                for (auto &atom: residue.atoms) {
                    atom_num++;
                }
            }
        }
        return atom_num;
    }


    void push(const Chain &);
    Chain &operator [](int);
    const Chain &operator [](int) const;

    /* IO function */
    void print();
    void write(string);
    friend ostream &operator <<(ostream &, const Model &);
    
    /* attributes */
    string name = "none";
    vector<Chain> chains;

};

} /// namespace jian

#endif // Model_H

