#pragma once

#include <iostream>
#include <set>
#include "pdb.hpp"
#include "geom.hpp"

namespace jian {

class Connect {
public:
    Model _model_1;
    Model _model_2;
    int _hinge_size = 2;
    std::set<std::string> _superposed_atoms{"C5*", "O3*", "C1*"};

    Model operator ()(const Model &a, const Model &b, int a1, int a2) {
        _model_1 = a;
        _model_2 = b;

        int len_1 = num_residues(_model_1);
        int len_2 = num_residues(_model_2);
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
        for (auto &chain : _model_1) {
            for (auto &residue : chain) {
                if (indices1.count(num_res)) {
                    for (auto &atom : residue) {
                        if (_superposed_atoms.count(atom.name)) {
                            for (int k = 0; k < 3; k++) x(temp, k) = atom[k];
                            temp++;
                        }
                    }
                }
                num_res++;
            }
        }
        num_res = 0;
        temp = 0;
        for (auto &chain : _model_2) {
            for (auto &residue : chain) {
                if (indices2.count(num_res)) {
                    for (auto &atom : residue) {
                        if (_superposed_atoms.count(atom.name)) {
                            for (int k = 0; k < 3; k++) y(temp, k) = atom[k];
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
        Model temp_rna;
        Chain temp_chain;
        num_res = 0;
        for (auto &chain : _model_1) {
            for (auto &residue : chain) {
                if (num_res <= a1 - _hinge_size || num_res >= a2 + _hinge_size) {
                    temp_chain.push_back(residue);
                } else if (num_res == a1) {
                    for (auto &chain2 : _model_2) {
                        for (auto &residue2 : chain2) {
                            temp_chain.push_back(residue2);
                        }
                    }
                }
                num_res++;
            }
        }
        temp_rna.push_back(temp_chain);
        return temp_rna;
    }

    template<typename T, typename L>
    void translate_model(T &model, const L &point) {
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) geom::translate(atom, point);
    }

    template<typename T, typename M>
    void rotate_model(T &model, const M &matrix) {
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) geom::rotate(atom, matrix);
    }

};

}

