#include "Transform.h"

namespace jian {

namespace nuc3d {

Model Transform::operator() (string type, string seq) {
    if (upper(type) == "RNA") {
        return to_rna(seq);
    } else if (upper(type) == "DNA") {
        return to_dna(seq);
    }
}

Model Transform::to_rna(string seq) {
    if (_model.res_nums() != seq.size()) 
        die("Transform::to_rna(string) error! Sequence's lenght isn't appropriate");

    int i = 0;
    for (auto &&chain: _model.chains) {
        for (auto &&residue: chain.residues) {
            cvt(residue, std::string() + seq[i]);
            i++;
        }
    }
    return _model;
}

Model Transform::to_dna(string seq) {
    if (_model.res_nums() != seq.size()) 
        die("Transform::to_rna(string) error! Sequence's lenght isn't appropriate");

    int i = 0;
    for (auto &&chain: _model.chains) {
        for (auto &&residue: chain.residues) {
            cvt(residue, std::string("D") + seq[i]);
            i++;
        }
    }
    return _model;
}

} /// namespace nuc3d

} /// namespace jian





