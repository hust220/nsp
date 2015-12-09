#ifndef JIAN_NUC3D_LM2_H
#define JIAN_NUC3D_LM2_H

#include <pdb/util.h>
#include <nuc2d/N2D.h>
#include <dg/DG.h>
#include <geom/move.h>
#include <geom/rotate.h>
#include <geom/util.h>
#include "Connect.h"
#include "Transform.h"
#include "JobInf.h"

namespace jian {
namespace nuc3d {

class LM2 : public virtual JobInf {
public:
    LM2(const Par &);
    Model operator ()();
    void read_pars();
    void set_native();
    void set_bound();
    void init();
    Model run();
    std::string parse_molecule_type(std::string seq);

    void set_helix_anchors();
    void set_helix_anchors(nuc2d::loop *src);
    void set_helix_coords(nuc2d::loop *src);
    Model to_all_atom(const MatrixXf &);
    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> get_fragments();
    MatrixXf make_helix_strand(MatrixXf, MatrixXf, int);
    std::pair<MatrixXf, std::deque<Residue>> make_helix(const MatrixXf &, int);
    std::pair<MatrixXf, std::deque<Residue>> make_frag(const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> &, MatrixXf, MatrixXf);

    template<typename List> std::deque<Residue> find_best_frag_model(List &&list) {
        typedef struct {
            bool operator ()(const std::pair<int, double> &p1, const std::pair<int, double> &p2) {return p1.second < p2.second;}
        } Compare;
        std::priority_queue<std::pair<int, double>, std::deque<std::pair<int, double>>, Compare> rank;
        int size = std::distance(std::begin(list), std::end(list));
        if (size == 3) {
            for (int i = 0; i < _frag_3_dists.size(); i++) 
                rank.push(std::make_pair(i, geometry::distance(_frag_3_dists[i], std::forward<List>(list))));
            std::cout << _frag_par_path + "/" + _frag_3_names[rank.top().first] + ".pdb" << std::endl;
            return get_residues_from_file(_frag_par_path + "/" + _frag_3_names[rank.top().first] + ".pdb");
        } else if (size == 10) {
            for (int i = 0; i < _frag_5_dists.size(); i++) 
                rank.push(std::make_pair(i, geometry::norm(_frag_5_dists[i], std::forward<List>(list), 10)));
            std::cout << _frag_par_path + "/" + _frag_5_names[rank.top().first] + ".pdb" << std::endl;
            return get_residues_from_file(_frag_par_path + "/" + _frag_5_names[rank.top().first] + ".pdb");
        } else {
            throw "JIAN::NUC3D::LM2::find_best_frag_model error!";
        }
    }


    DG _dg;
    std::map<int, int> _m1, _m2;
    std::string _job;
    MatrixXf _native_helices;
    MatrixXf _helices;
    MatrixXf _native_scaffold;
    MatrixXf _scaffold;
    Model _native_model;

    map<string, MatrixXf> _mono_nuc_pars;
    map<string, MatrixXf> _adj_nuc_pars;
    map<string, MatrixXf> _aa_pars;

//    void set_base_pairs(nuc2d::loop *);

    MatrixXf get_helix_par(const R5P &);
    R5P get_helix(const nuc2d::helix &);
    R5P create_helix(const string &);

    nuc2d::N2D _ss_tree;
    std::vector<std::vector<nuc2d::res>> _helix_anchors;
    std::map<nuc2d::loop *, MatrixXf> _helix_coords;
    std::map<nuc2d::loop *, MatrixXf> _nuc_coords;
    std::map<nuc2d::loop *, MatrixXf> _atom_coords;
    MatrixXf _bound;
    DG dg;
    int _atom_nums_per_nuc = 5;
    double _err_radius = 0.001;
    Model _helix_file;
    int _view = 0;
    std::string _frag_par_path;
    std::vector<std::string> _frag_3_names;
    std::vector<std::string> _frag_4_names;
    std::vector<std::string> _frag_5_names;
    std::vector<std::array<double,  3>> _frag_3_dists;
    std::vector<std::array<double,  6>> _frag_4_dists;
    std::vector<std::array<double, 10>> _frag_5_dists;

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

