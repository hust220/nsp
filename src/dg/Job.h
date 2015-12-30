#ifndef JIAN_DG_JOB_H
#define JIAN_DG_JOB_H

#include "../util/std.h"
#include "../util/mat.h"
#include "../util/Log.h"
#include "DihBound.h"

namespace jian {
namespace dg {

class Job {
public:
    using DihEntry = std::vector<int>;
    using DistType = MatrixXf;
    using DihType = DihBound<DihEntry>;

    std::mt19937 _rand_engine{11};
    std::uniform_real_distribution<double> _unif_distr{0, 1};

    int len = 0;

    MatrixXf bound;
    double _min_dist = 0;
    DihBound<DihEntry> _dih_bound;
    MatrixXf d;
    MatrixXf m;
    MatrixXf c;
    MatrixXf g;
    double g2 = 0;
    double _dist_en = 0;
    double _dih_en = 0;

    Log log;

    void seed(int n) {
        _rand_engine.seed(n);
    }

    double rand() {
        return _unif_distr(_rand_engine);
    }
};

} // namespace dg
} // namespace jian

#endif

