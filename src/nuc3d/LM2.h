#ifndef JIAN_NUC3D_LM2_H
#define JIAN_NUC3D_LM2_H

#include <pdb/util.h>
#include <nuc2d/util.h>
#include <dg/DG.h>
#include "Connect.h"

namespace jian {
namespace nuc3d {

class LM2 {
public:
    LM2(string type = "RNA");
    std::vector<Model> operator ()(string seq, string ss, string constraint_file = "", int num = 1);

    void set_helix_anchors(nuc2d::loop *src);
    void set_helix_coords(nuc2d::loop *src);
    Model to_all_atom(const MatrixXf &);
    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> get_fragments(const std::vector<int> &);
    MatrixXf get_frag_coords(const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> &, MatrixXf, MatrixXf);
    MatrixXf make_helix_strand(MatrixXf, MatrixXf, int);
    MatrixXf make_helix(const MatrixXf &, int);


    map<string, MatrixXf> _mono_nuc_pars;
    map<string, MatrixXf> _adj_nuc_pars;
    map<string, MatrixXf> _aa_pars;

    void init();
//    void set_base_pairs(nuc2d::loop *);

    MatrixXf get_helix_par(const R5P &);
    R5P get_helix(const nuc2d::helix &);
    R5P create_helix(const string &);

    nuc2d::N2D _ss_tree;
    MatrixXf _scaffold;
    std::vector<std::vector<nuc2d::res>> _helix_anchors;
    std::map<nuc2d::loop *, MatrixXf> _helix_coords;
    std::map<nuc2d::loop *, MatrixXf> _nuc_coords;
    std::map<nuc2d::loop *, MatrixXf> _atom_coords;
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

