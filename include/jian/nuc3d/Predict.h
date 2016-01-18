#ifndef JIAN_NUC3D_PREDICT_H
#define JIAN_NUC3D_PREDICT_H

#include "../nuc2d/N2D.h"
#include "../nuc2d/util.h"
#include "../pdb/util.h"
#include "Assemble.h"
#include "triass/TriAss.h"

namespace jian {
namespace nuc3d {

class Predict : public virtual JobInf {
public:
    Par _par;

    Predict(const Par &par) : JobInf(par), _par(par) {}

    void operator ()() {
        predict();
    }

    void predict() {
        boost::timer t;
        log.clear();
        log("=========================================================\n",
            "New Job: ", _name, '\n',
            "Time: ", Time::time(), '\n',
            "Sequence: ", _seq, '\n',
            "2D Structure: ", _ss, '\n',
            "Number: ", _num, '\n',
            "Method: ", _method, '\n',
            "----------------------------------------\n");

        if (_method == "assemble") {
            Assemble assemble(_par);
            assemble();
        } else if (_method == "tri-assemble") {
            triass::TriAss<DNA> tri_ass(_par);
            tri_ass();
        }

        log("----------------------------------------\n",
            "Finish Time: ", Time::time(), '\n',
            "Time elapsed: ", t.elapsed(), "s\n",
            "=========================================================\n\n");
    }

};

} // namespace nuc3d
} // namespace jian

#endif

