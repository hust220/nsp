#ifndef JIAN_SCORING_TRAIN_H
#define JIAN_SCORING_TRAIN_H

#include "DistAnal.h"
#include "DihAnal.h"

namespace jian {
namespace scoring {

class Train {
public:
    Train(Par);
    ~Train();
    void operator ()();
//    void dih(char *);
    void print_par(const MatrixXd &);
    
    DistAnal *_dist_anal;
    std::string _file_list;
    int _cutoff = 20;
    double _bin = 0.3;
    int _bins = 67;
    string _par_dist;
    string _par_dih;
};

} /// namespace scoring
} /// namespace jian

#endif

