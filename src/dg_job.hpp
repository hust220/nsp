#pragma once

#include "rand.hpp"
#include "log.hpp"
#include "matrix.hpp"
#include "dg_dih_bound.hpp"

BEGIN_JN

namespace dg {

struct Job
{
    using DihEntry = Vi;
    using DistBoundType = Mat;
    using DihBoundType = DihBound;

    int len = 0;

    Log log;

    Mat bound; /**< distance bound matrix. */
    DihBound _dih_bound; /**< dihedral bound matrix. */

    double _min_dist = 5;
    Mat d, m, c, g;
    double g2 = 0;
    double _dist_en = 0;
    double _dih_en = 0;
    int _num_mc_steps = 2000;

};

} // namespace dg

END_JN

