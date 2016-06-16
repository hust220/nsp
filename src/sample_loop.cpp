#include "nsp.hpp"
#include <jian/nuc3d/SplitModel.hpp>
#include <jian/nuc3d/SampleLoop.hpp>

namespace jian {

REGISTER_NSP_COMPONENT(sample_loop) {
    loop *l;
    std::shared_ptr<Model> model;
    SplitModel split_model(par["model"][0], par["seq"][0], par["ss"][0]);
    LOOP_TRAVERSE(split_model._ss_tree.head(),
        if (L->num_sons() == 3) {
            l = L;
            model = std::make_shared<Model>(split_model._loops[l]);
        }
    );
    SampleLoop sample_loop(l, *model);
    for (int i = 0; i < JN_INT(par["num"][0]); i++) {
        write_pdb(sample_loop(), par["out"][0] + JN_STR(i+1) + ".pdb");
    }
}

} // namespace jian

