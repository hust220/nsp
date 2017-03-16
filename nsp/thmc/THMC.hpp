#pragma once

#include <memory>
#include <Eigen/Dense>
#include <map>
#include <deque>
#include <numeric>
#include <vector>
#include "TModule.hpp"
#include "../pdb.hpp"
#include "jian/geom.hpp"
#include "jian/pp.hpp"
#include "../nuc3d/TSP.hpp" 
#include "../mcsm.hpp"
#include "jian/utils/Factory.hpp"
#include "jian/utils/Par.hpp"
#include "jian/utils/rand.hpp"
#include "jian/utils/Env.hpp"
#include "jian/utils/file.hpp"

BEGIN_JN
namespace nuc3d {
namespace triple {

using fac_t = Factory<TModule::cons_t>;

class THMC : public MCSM {
public:
    using Res = struct {char seq; char ss; int num;};
    using RelatedResidues = Vector<SP<Deque<int>>>;

    Tree _tree;
    Deque<TModule *> d_modules;
    int d_mc_selected_index;
    RelatedResidues d_mc_related_residues;
    RelatedResidues d_mc_unrelated_residues;

    THMC() = default;

    void init(const Par &par) {
        MCSM::init(par);
    }

    ~THMC();

    void print_related_residues(const RelatedResidues &r);

    void print_modules();

    void set_modules();

    void ss_to_tree();

    void build_initial_scaffold();

    Chain build_helix(int len);

    Chain load_triple_helix(int n);

    void set_coords_residue(Mat &c1, int m, const Residue &r);

    Chain connect_triple_helix(Chain &c1, Chain &c2);

    void shrink_to_fit(const Chain &c);

    void print_tuple(const Tuple &tuple);

    void print_helix(const Tuples &helix);

    void print_tree();

    // mc-related functions

    void init_for_mc();

    virtual void before_run();

    virtual void mc_sample();

    virtual void mc_select();

    virtual bool is_selected(const int &i) const;

    virtual Vec rotating_center() const;
};

} // namespace triple
} // namespace nuc3d
END_JN


