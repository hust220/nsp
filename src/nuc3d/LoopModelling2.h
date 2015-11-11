/**
 * @file LoopModelling22.cpp
 * @brief Modelling loop by using distance geometry algorithm.
 *
 * Each nucleotide is represented by 5 atoms.
 *   A: C5* O3* C1* N6 C2
 *   U: C5* O3* C1* O2 O4
 *   G: C5* O3* C1* O6 N2
 *   C: C5* O3* C1* O2 N4
 *
 * @author wj_hust08@163.com
 * @version 1.0
 * @date 2015-5-25
*/

#ifndef LOOPMODELLING2_H
#define LOOPMODELLING2_H

#include <pdb/util.h>
#include <nuc2d/util.h>
#include <dg/DG.h>
#include "Connect.h"

namespace jian {

namespace nuc3d {

class LoopModelling2 {
public:
    LoopModelling2(string type = "RNA");
    Model operator ()(string seq, string ss, string constraint_file = "", int num = 1);
    Model to_all_atom(const MatrixXf &);

    map<string, MatrixXf> _mono_nuc_pars;
    map<string, MatrixXf> _adj_nuc_pars;
    map<string, MatrixXf> _aa_pars;

    void init();
    void set_base_pairs(nuc2d::loop *);

    MatrixXf get_helix_par(const R5P &);
    R5P get_helix(const nuc2d::helix &);
    R5P create_helix(const string &);

    string _seq;
    string _ss;
    MatrixXf _bound;
    DG dg;
    int _atom_nums_per_nuc = 5;
    double _err_radius = 0.001;
    string _lib;
    string _type = "RNA";
    int _view = 0;

    /////////////////////////////////////
    /// Constraints
    void set_constraints();
    MatrixXf apply_constraints(const MatrixXf &);
    std::vector<std::vector<std::vector<int>>> read_constraints_file();
    MatrixXf read_constraints_pos(std::string, int);
    MatrixXf read_constraints_par(std::string, int);

    string _constraint_file;
    //////////////////////////////////////

};

} // namespace nuc3d

} // namespace jian

#endif

