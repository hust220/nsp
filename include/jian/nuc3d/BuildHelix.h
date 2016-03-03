#ifndef JIAN_NUC3D_BUILDHELIX
#define JIAN_NUC3D_BUILDHELIX

#include <nuc2d/N2D.h>
#include "Connect.h"

namespace jian {
namespace nuc3d {

class BuildHelix {
public:
    std::string _lib;
    std::string _type = "RNA";

    BuildHelix() {
        _lib = env("NSP") + "/" + _type;
    }

    Model operator ()(nuc2d::helix *s) {
        std::string seq1, seq2;
        for (nuc2d::bp *b = s->head; b != NULL; b = b->next) {
            seq1 += b->res1.name;
            seq2 += b->res2.name;
        }
        return (*this)(seq1 + seq2);
    }

    Model operator ()(std::string seq) {
        if (seq.size() < 2 || seq.size() % 2 == 1) {
            die("jian::nuc3d::BuildHelix::operator (std::string) error! Unreasonable length.\nSequence: " + seq);
        } else if (seq.size() == 2) {
            std::string file_name = _lib + "/basepair/" + seq + ".pdb";
            std::ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XX.pdb";
            }
            ifile.close();
            return Model(file_name);
        } else if (seq.size() == 4) {
            std::string file_name = _lib + "/basepair/" + seq + ".pdb";
            std::ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XXXX.pdb";
            }
            ifile.close();
            return Model(file_name);
        } else {
            std::string file_name = _lib + "/basepair/" + seq.substr(0, 2) + seq.substr(seq.size() - 2, 2) + ".pdb";
            std::ifstream ifile(file_name.c_str());
            if (!ifile) {
                file_name = _lib + "/basepair/XXXX.pdb";
            }
            ifile.close();
            Connect connect;
            connect._hinge_size = 1;
            return connect(Model(file_name), (*this)(seq.substr(1, seq.size() - 2)), 1, 2);
        }
    }

};

} /// namespace nuc3d
} /// namespace jian

#endif

