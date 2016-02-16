#ifndef JIAN_DG_JOB_H
#define JIAN_DG_JOB_H

#include "../etl.h"
#include "DihBound.h"

namespace jian {
namespace dg {

class Job : public virtual Rand {
public:
    using Mat = MatrixXd;
    using DihEntry = std::vector<int>;
    using DistBoundType = Mat;
    using DihBoundType = DihBound<DihEntry>;

    int len = 0;

    Mat bound;
    double _min_dist = 0;
    DihBound<DihEntry> _dih_bound;
    Mat d, m, c, g;
    double g2 = 0;
    double _dist_en = 0;
    double _dih_en = 0;

    Log log;

};

} // namespace dg
} // namespace jian

#endif

