#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "../pdb/util.h"
#include "../util/Log.h"

namespace jian {
namespace nuc3d {

class Transform {
public:
    Model _model;
    Convert cvt;
    Log log;

    Transform() {}

    Transform(const Model &model) {
        _model = model;
    }

    Transform(Model &&model) {
        std::swap(_model, model);
    }

    Model &model() {
        return _model;
    }

    const Model &model() const {
        return _model;
    }

    void model(const Model &model) {
        _model = model;
    }

    Model operator() (string type, string seq) {
        if (upper(type) == "RNA") {
            return to_rna(seq);
        } else if (upper(type) == "DNA") {
            return to_dna(seq);
        }
    }

    Model to_rna(string seq) {
        if (_model.res_nums() != seq.size()) {
            for (auto &&chain: _model.chains) for (auto &&res: chain.residues) std::cout << res.name;
            std::cout << std::endl;
            die("Transform::to_rna(string) error! Sequence's length isn't appropriate.\nNumbers of residues: " + boost::lexical_cast<std::string>(_model.res_nums()) + ".\nLength of sequence: " + boost::lexical_cast<std::string>(seq.size()) + ".\nSequence: " + seq + ".");
        }
        int i = 0;
        for (auto &&chain: _model.chains) for (auto &&residue: chain.residues) {
            cvt(residue, std::string() + seq[i]);
            i++;
        }
        return _model;
    }

    Model to_dna(string seq) {
        if (_model.res_nums() != seq.size()) die("Transform::to_rna(string) error! Sequence's lenght isn't appropriate");

        int i = 0;
        for (auto &&chain: _model.chains) for (auto &&residue: chain.residues) {
            cvt(residue, std::string("D") + seq[i]);
            i++;
        }
        return _model;
    }

};


}


}





#endif 

