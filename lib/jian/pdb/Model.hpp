#pragma once

#include "Chain.hpp"
#include <string>
#include "../utils/traits.hpp"

#define EACH_INDEX_RES_HELPER(n, m, c) ({\
    int N_CHAIN##n = 0; \
    int N_RES##n = 0; \
    for (auto && CHAIN##n : m) {\
        for (auto && RES##n : CHAIN##n) {\
            c;\
            N_RES##n++;\
        }\
        N_CHAIN##n++;\
    }\
    N_RES##n;\
})
#define EACH_RES1(m, c) EACH_INDEX_RES_HELPER(1, m, c)
#define EACH_RES2(m, c) EACH_INDEX_RES_HELPER(2, m, c)
#define EACH_RES3(m, c) EACH_INDEX_RES_HELPER(3, m, c)
#define EACH_RES4(m, c) EACH_INDEX_RES_HELPER(4, m, c)
#define EACH_RES(m, c)  EACH_INDEX_RES_HELPER( , m, c)

#define EACH_ATOM_HELPER(n, m, c) ({\
    int N_ATOM##n = 0;\
    int N_RES##n = 0; \
    int N_CHAIN##n = 0; \
    for (auto && CHAIN##n : m) {\
        for (auto && RES##n : CHAIN##n) {\
            for (auto && ATOM##n : RES##n) {\
                c;\
                N_ATOM##n++;\
            }\
            N_RES##n++;\
        }\
        N_CHAIN##n++;\
    }\
    N_ATOM##n;\
})
#define EACH_ATOM1(m, c) EACH_ATOM_HELPER(1, m, c)
#define EACH_ATOM2(m, c) EACH_ATOM_HELPER(2, m, c)
#define EACH_ATOM3(m, c) EACH_ATOM_HELPER(3, m, c)
#define EACH_ATOM4(m, c) EACH_ATOM_HELPER(4, m, c)
#define EACH_ATOM(m, c)  EACH_ATOM_HELPER( , m, c)

namespace jian {

class Model : public std::deque<Chain> {
public:
    std::string name = "unknown";
    std::string type = "unknown";
    int num = 1;
};

std::string seq(const Model &model);
int num_residues(const Model &model);
int num_atoms(const Model &model);
bool is_empty(const Model &model);

template<typename T> 
auto sub(const Model &model, T &&t) {
    Model m;
    int res_num = 0; 
    for (auto &&chain: model) {
        Chain temp_chain; 
        temp_chain.name = chain.name;
        for (auto &&res: chain) {
            if (std::find(std::begin(t), std::end(t), res_num) != std::end(t)) temp_chain.push_back(res);
            res_num++;
        }
        if (!temp_chain.empty()) m.push_back(temp_chain);
    }
    return m;
}

template<template<typename...> class L = std::deque>
auto residues(const Model &model) {
    L<Residue> v;
    for (auto &&chain: model) for (auto &&res: chain) v.push_back(res);
    return v;
}

template<typename T, template<typename...> class L = std::deque> 
auto residues(const Model &model, T &&ls) {
    L<Residue> v;
    int i = 0; for (auto &&chain: model) for (auto &&res: chain) {
        if (std::count(std::begin(ls), std::end(ls), i)) v.push_back(res);
        i++;
    }
    return v;
}

template<typename T>
uniform_const_t<Residue, T> &residue(T &&model, int n) {
    int i = 0; for (auto &&chain: model) for (auto &&res: chain) {
        if (i == n) return res;
        i++;
    }
    throw "JIAN::MODEL::residue(int) error! Residue index out of range.";
}

template<typename T>
auto residues_to_model(T &&ls) {
    Model model; model.resize(1);
    for (auto &&res : ls) model[0].push_back(res);
    return model;
}

} // namespace jian

