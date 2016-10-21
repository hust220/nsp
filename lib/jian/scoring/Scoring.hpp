#pragma once

#include "DistAnal.hpp"
#include "DihAnal.hpp"

namespace jian {

//class DistAnal;
//class DihAnal;

class Scoring {
public:
    DistAnal * m_dist_anal = NULL;
    DihAnal  * m_dih_anal = NULL;

    std::string m_par_dist;
    std::string m_par_dih;

    double m_bin_dist = 0.3;
    double m_bin_dih  = 4.5;
    double m_cutoff   = 20;

    double m_score_dist = 0;
    double m_score_dih = 0;
    double m_score = 0;

    double m_constant = 27.1118;
    double m_weight_dist = 0.433513;
    double m_weight_dih = 1.59348;


    void init();
    ~Scoring();
    Scoring & run(const Chain &);
	Scoring & train(const Chain &);
};

} // namespace jian

