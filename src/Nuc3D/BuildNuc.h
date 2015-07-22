#ifndef BUILD_NUC_H
#define BUILD_NUC_H

#include "../DG.h"
#include "../Pdb.h"

namespace jian {

namespace nuc3d {

class BuildNuc {
public:
    BuildNuc(string type = "RNA");
    Residue operator() (const string &name, const MatrixXf &scaffold);
    Residue make_residue(const string &name, const MatrixXf &coords);

    int _view = 0;

private:
    string _type = "RNA";
    string _lib;
    map<string, MatrixXf> _base_aa_par;
    map<string, MatrixXf> _base_cg_par;
    map<string, MatrixXf> _phos_sugar_par;
    map<string, vector<string>> _atom_names;
};

} /// namespace nuc3d

} /// namespace jian





#endif

