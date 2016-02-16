#ifndef CONNECT_H
#define CONNECT_H

#include "../pdb.h"

namespace jian {
namespace nuc3d {

class Connect {
public:
    using Mat = MatrixXd;

    Model _model_1;
    Model _model_2;
    int _hinge_size = 2;
    set<string> _superposed_atoms{"C5*", "O3*", "C1*"};

    Model operator ()(const Model &a, const Model &b, int a1, int a2) {
        _model_1 = a;
        _model_2 = b;

        int len_1 = _model_1.res_nums();
        int len_2 = _model_2.res_nums();
        int superposed_atom_nums = _superposed_atoms.size() * _hinge_size * 2;
        std::set<int> indices1, indices2;
        for (int i = 0; i < _hinge_size; i++) {
            indices1.insert(a1 - i);
            indices1.insert(a2 + i);
            indices2.insert(0 + i);
            indices2.insert(len_2 - 1 - i);
        }

        Mat x(superposed_atom_nums, 3), y(superposed_atom_nums, 3);
        int num_res = 0;
        int temp = 0;
        for (auto &chain : _model_1.chains) {
            for (auto &residue : chain.residues) {
                if (indices1.count(num_res)) {
                    for (auto &atom : residue.atoms) {
                        if (_superposed_atoms.count(atom.name)) {
                            x(temp, 0) = atom.x;
                            x(temp, 1) = atom.y;
                            x(temp, 2) = atom.z;
                            temp++;
                        }
                    }
                }
                num_res++;
            }
        }
        num_res = 0;
        temp = 0;
        for (auto &chain : _model_2.chains) {
            for (auto &residue : chain.residues) {
                if (indices2.count(num_res)) {
                    for (auto &atom : residue.atoms) {
                        if (_superposed_atoms.count(atom.name)) {
                            y(temp, 0) = atom.x;
                            y(temp, 1) = atom.y;
                            y(temp, 2) = atom.z;
                            temp++;
                        }
                    }
                }
                num_res++;
            }
        }

        auto sp = geom::suppos(x, y);

        translate_model(_model_1, -sp.c1);
        translate_model(_model_2, -sp.c2);
        rotate_model(_model_1, sp.rot);

        // construct a new RNA
        RNA temp_rna;
        Chain temp_chain;
        num_res = 0;
        for (auto &chain : _model_1.chains) {
            for (auto &residue : chain.residues) {
                if (num_res <= a1 - _hinge_size || num_res >= a2 + _hinge_size) {
                    temp_chain.residues.push_back(residue);
                } else if (num_res == a1) {
                    for (auto &chain2 : _model_2.chains) {
                        for (auto &residue2 : chain2.residues) {
                            temp_chain.residues.push_back(residue2);
                        }
                    }
                }
                num_res++;
            }
        }
        temp_rna.chains.push_back(temp_chain);
        return temp_rna;
    }

    template<typename T, typename L>
    void translate_model(T &model, const L &point) {
        for (auto &chain : model.chains) {
            for (auto &residue : chain.residues) {
                for (auto &atom : residue.atoms) {
                    geom::translate(atom, point);
                }
            }
        }
    }

    template<typename T, typename M>
    void rotate_model(T &model, const M &matrix) {
        for (auto &chain : model.chains) {
            for (auto &residue : chain.residues) {
                for (auto &atom : residue.atoms) {
                    geom::rotate(atom, matrix);
                }
            }
        }
    }

};

} // namespace nuc3d
} // namespace jian

#endif //CONNECT_H







