#ifndef SCORE_H
#define SCORE_H

#include "DistAnal.h"
#include "DihAnal.h"

namespace jian {
namespace scoring {

class Score {
public:
//    Score();
//    Score(char *);
    Score(Par par);
    ~Score();
    void operator ()();

    DistAnal * _dist_anal;
    DihAnal * _dih_anal;

    std::string _par_dist_obs;
    std::string _par_dist_ref;
    std::string _par_dih;
    std::string _ref_state = "average";
    std::string _lib;
    std::string _file_list;
    std::string _input;

    double _dist_bin = 0.3;
    double _dih_bin = 4.5;
    int _cutoff = 20;

    double _constant = 27.1118;
    double _dist_weight = 0.433513;
    double _dih_weight = 1.59348;

    std::string _type = "RNA";

};

} /// namespace scoring
} /// namespace jian

#endif

