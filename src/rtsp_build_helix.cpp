#include <iostream>
#include <fstream>
#include <string>
#include "rtsp_build_helix.hpp"
#include "rtsp_mutate.hpp"
#include "env.hpp"
#include "geom.hpp"

BEGIN_JN

static void res_set_mat(Mat &m, const Residue &r, Int n = 0) {
    static Vector<Str> atoms {"C5*", "O3*", "C1*"};
    Array<Bool, 3> flag {false, false, false};
    for (auto && atom : r) {
        auto it = std::find(atoms.begin(), atoms.end(), atom.name);
        if (it != atoms.end()) {
            Int d = std::distance(atoms.begin(), it);
            flag[d] = true;
            for (Int i = 0; i < 3; i++) {
                m(n+d, i) = atom[i];
            }
        }
    }
    for (Int i = 0; i < 3; i++) if (!flag[i]) throw to_str("build_helix.cpp line: ", __LINE__, ", atom: ", atoms[i], " not found!");
}

static auto bp_sp(const Residue &r1, const Residue &r2, const Residue &s1, const Residue s2) {
    Mat m1(6, 3), m2(6, 3);
    res_set_mat(m1, r1);
    res_set_mat(m1, r2, 3);
    res_set_mat(m2, s1);
    res_set_mat(m2, s2, 3);
    return geom::suppos(m1, m2);
}

static Chain connect(Chain &&c1, Chain &&c2) {
    Int l1 = size(c1), l2 = size(c2);
    if (l1 < 4 || l2 < 2) throw to_str("build_helix.cpp line: ", __LINE__, ", l1: ", l1, ", l2: ", l2);
    Int n = l1/2;
    auto sp = bp_sp(c1[n-1], c1[n], c2.front(), c2.back());
    for (auto && res : c1) for (auto && atom : res) sp.apply(atom);
    Chain c(std::move(c2));
    c.push_front(std::move(c1.front()));
    c.push_back(std::move(c1.back()));
    return c;
}

static Str chain_seq(const Chain &c) {
    Str seq;
    for (auto && res : c) {
        seq += res.name.back();
    }
    return seq;
}

Chain build_helix(S seq) {
    S lib = Env::lib() + "/RNA";
    if (seq.size() < 2 || seq.size() % 2 == 1) {
        throw std::string("jian::BuildHelix::operator (std::string) error! Unreasonable length.\nSequence: ") + seq;
    } else if (seq.size() == 2) {
        S file_name = lib + "/basepair/" + seq + ".pdb";
        std::ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = lib + "/basepair/XX.pdb";
        }
        ifile.close();
        return mol_read_to<Chain>(file_name);
    } else if (seq.size() == 4) {
        S file_name = lib + "/basepair/" + seq + ".pdb";
        std::ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = lib + "/basepair/XXXX.pdb";
        }
        ifile.close();
        return mol_read_to<Chain>(file_name);
    } else {
        S file_name = lib + "/basepair/" + seq.substr(0, 2) + seq.substr(seq.size() - 2, 2) + ".pdb";
        std::ifstream ifile(file_name.c_str());
        if (!ifile) {
            file_name = lib + "/basepair/XXXX.pdb";
        }
        ifile.close();
        Chain && c = connect(mol_read_to<Chain>(file_name), build_helix(seq.substr(1, seq.size() - 2)));
        if (chain_seq(c) != seq) {
            return mutate(c, seq, "RNA");
        }
        else {
            return c;
        }
    }
}

END_JN

