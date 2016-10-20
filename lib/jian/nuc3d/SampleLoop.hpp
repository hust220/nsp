#pragma once

#include <iostream>
#include <deque>
#include <string>
#include "../utils/rand.hpp"
#include "../pdb/Model.hpp"
#include "../nuc2d/SSTree.hpp"
#include "../geom.hpp"

namespace jian {

class SampleLoop {
public:
    using Ref = struct {int beg, end, len;};
    using RandPars = struct {int i_trans; double v_trans; int i_rot; double v_rot;};

    std::deque<std::array<double, 3>> _centers;
    std::deque<Ref> _refs;
    int _num_helices;
    std::string _seq;
    std::string _ss;
    Model _model;

    void set_refs() {
        SSTree ss_tree; ss_tree.make(_seq, _ss);
        LOOP_TRAVERSE(ss_tree.head(),
            if (L->has_helix()) {
                _refs.push_back({L->s.head->res1.num-1, L->s.head->res2.num-1, L->s.len()});
            }
        );
    }

    void set_centers() {
        for (int i = 0; i < _num_helices; i++) {
            _centers.push_back({0, 0, 0});
        }
        EACH_RES(_model,
            for (int i = 0; i < _num_helices; i++) {
                if ((N_RES >= _refs[i].beg && N_RES < _refs[i].beg + _refs[i].len) ||
                    (N_RES > _refs[i].end - _refs[i].len && N_RES <= _refs[i].end)) {
                    auto &&atom = RES["C4*"];
                    for (int j = 0; j < 3; j++) {
                        _centers[i][j] += atom[j];
                    }
                }
            }
        );
        for (int i = 0; i < _num_helices; i++) {
            for (int j = 0; j < 3; j++) {
                _centers[i][j] /= _refs[i].len * 2;
            }
        }
    }

    SampleLoop(loop *l, const Model &model) {
        _seq = l->seq(); 
        _ss = NucSS::lower_ss(l->ss());
        _model = model;
        _num_helices = l->num_branches();
        set_refs();
        set_centers();
    }

    RandPars rand_pars() {
        return RandPars {
            int(rand() * 3), (rand() - 0.5) * 2,
            int(rand() * 3), (rand() - 0.5) * PI / 6
        };
    }

    template<typename T>
    void apply_pars(int index, T &&pars) {
        int &beg = _refs[index].beg;
        int &end = _refs[index].end;
        int &len = _refs[index].len;
        auto &&mat = geom::rot_mat(pars.i_rot, pars.v_rot);
        EACH_RES(_model,
            if ((N_RES >= beg && N_RES < beg + len) || (N_RES > end - len && N_RES <= end)) {
                EACH(atom, RES,
                    FOR((i, 3), atom[i] -= _centers[index][i]);
                    geom::rotate(atom, mat);
                    _centers[index][pars.i_trans] += pars.v_trans;
                    FOR((i, 3), atom[i] += _centers[index][i]);
                );
            }
        );
    }

    Model operator ()() {
        int index = int(rand() * _num_helices);
        auto pars = rand_pars();
        apply_pars(index, pars);
        return _model;
    }

};

} // namespace jian

