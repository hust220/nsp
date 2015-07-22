#ifndef SCORE_H
#define SCORE_H

#include "DistAnal.h"
#include "DihAnal.h"

namespace jian {

class Score {
public:
    Score();
    Score(char *);
    double operator ()(const RNA &);

    Obj<DistAnal> distAnal;
    Obj<DihAnal> dihAnal;

    string par_dist_obs;
    string par_dist_ref;
    string dihPar;
    string reference_state = "average";

    double _dist_bin = 0.3;
    double _dih_bin = 4.5;

    double _constant = 27.1118;
    double _distWeight = 0.433513;
    double _dihWeight = 1.59348;

    int cutoff = 20;
};

} /// namespace jian

#endif

