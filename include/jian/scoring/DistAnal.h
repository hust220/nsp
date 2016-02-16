#ifndef DISTANAL_H
#define DISTANAL_H

#include "Analysis.h"

namespace jian {
namespace scoring {

class DistAnal {
public:
    DistAnal(int, double, string = "average");

    void read_mol(const Model &);
    void train();
    double operator ()(const Model &);

    void read_obs_parm(string = "");
    void read_ref_parm(string = "");
    void init_obs_prob();
    void init_ref_prob();

    double score();

    VectorXf _scores = VectorXf::Zero(5);
    VectorXf _nuc_score;
    VectorXd _nuc_len;
    VectorXf _stacking_score;
    VectorXd _stacking_len;
    MatrixXf _pairing_score;
    MatrixXd _pairing_len;

    string _ref_state = "average";
    int _cutoff = 20;
    double _interval = 0.3;
    int _bins = 67;
    double _penalty = 0;

    VectorXd _num;
    VectorXd _type;
    VectorXd _ntLen;
    std::vector<std::vector<Point>> _list;
    int _len;

    /* 1: i-i+1; 2: i-i+2; 3: i-i+3; 4: i-i+n */
    MatrixXd _obs_parm;
    MatrixXd _ref_parm;
    MatrixXf _obs_prob;
    MatrixXf _ref_prob;
};

} /// namespace scoring
} /// namespace jian

#endif

