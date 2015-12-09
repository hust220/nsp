#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <nuc2d/N2D.h>
#include <nuc2d/util.h>
#include <pdb/util.h>
#include "FindTemplates.h"
#include "AssembleTemplates.h"

namespace jian {
namespace nuc3d {

class Assemble : public virtual AssembleTemplates {
public:
    Log log;

    Assemble(const Par &pars) : JobInf(pars) {
    }

    void operator ()() {
        assemble();
    }

    void assemble() {
        log("=========================================================");
        log("New Job: " + _name);
        log("Time: " + Time::time());
        boost::timer t;
        log("Sequence: " + _seq);
        log("2D Structure: " + _ss);
        log("Number: " + std::to_string(_num));
        log("----------------------------------------");

        log("Step 1: Construct 2D structure tree.");
        _n2d.hinge_base_pair_num = _hinge;
        _n2d(_seq, _ss);

        log("Step 2: Searching templates.");
        find_templates();

        log("Step 3: Assemble templates.");
        for (int i = 0; i < _num; i++) {
            log("assemble model " + std::to_string(i + 1) + "...");
            assemble_templates().write(_name + "-" + std::to_string(i) + ".pdb");
        }

        log("----------------------------------------");
        log("Finish Time: " + Time::time());
        log("Time elapsed: " + std::to_string(t.elapsed()) + "s");
        log("=========================================================");
        log("\n\n");
    }

};

} // namespace nuc3d
} // namespace jian

#endif

