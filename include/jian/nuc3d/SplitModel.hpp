#pragma once

#include "../pdb.hpp"
#include "../nuc2d.hpp"

namespace jian {

class SplitModel {
public:
    SSTree _ss_tree;
    std::map<loop *, Model> _helices;
    std::map<loop *, Model> _loops;
    std::shared_ptr<Model> _model;

    SplitModel(const std::string &name, const std::string &seq, const std::string &ss) {
        _model = std::make_shared<Model>(name);
        _ss_tree.make(seq, ss);
        LOOP_TRAVERSE(_ss_tree.head,
            parse_helix(L).parse_loop(L);
        );
    }

    SplitModel &parse_helix(loop *l) {
        if (l->has_helix()) {
            int len = l->s.len();
            std::vector<int> nums(len * 2);
            HELIX_EACH(l->s,
                nums[N_BP] = BP->res1.num - 1;
                nums[2 * len - 1 - N_BP] = BP->res2.num - 1;
            );
            _helices[l] = sub(*_model, nums);
        }
        return *this;
    }

    SplitModel &parse_loop(loop *l) {
        if (l->has_loop()) {
            std::deque<int> all_nums;
            LOOP_EACH(l, all_nums.push_back(RES->num - 1));
            _loops[l] = sub(*_model, all_nums);
        }
        return *this;
    }
};

} // namespace jian


