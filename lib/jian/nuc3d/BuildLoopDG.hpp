#pragma once

#include "../utils/Debug.hpp"
#include "../pdb.hpp"
#include "../nuc2d.hpp"
#include "../dg.hpp"
#include "HelixPar.hpp"

namespace jian {

Chain *build_chain_dg(std::string seq, std::string ss);

class BuildLoopDG {
public:
    Eigen::MatrixXd _dist_bound;
    DihBound _dih_bound;

    DG dg;
	HelixPar helix_par;

    void init(const std::string &seq, const std::string &ss, const std::string &c = "");

    Chain operator ()();

    void set_bound_loop(Eigen::MatrixXd &b, DihBound &d, loop *l);

    void set_bound_helix(Eigen::MatrixXd &b, DihBound &d, const helix &h);
};

} // namespace jian

