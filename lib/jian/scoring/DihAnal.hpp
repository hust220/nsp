#pragma once

#include <string>
#include <array>

BEGIN_JN

class Chain;

class DihAnal {
public:
    using Point = std::array<double, 3>;

    DihAnal(double = 1.0);
    ~DihAnal();

    void read_mol(const Chain &);
    void train();
    DihAnal &run(const Chain &);
    void read_parm(std::string);
    void initProb();
    void printParm();
    void printProb();
    double getScore();

    void delPoints();
    void initPoints(int);

    double interval;
    int bins;

    int *obsParm;
    double *obsProb;
    int *refParm;
    double *refProb;
    int len;
    Point **p, **o5_, **c5_, **c4_, **c3_, **o3_, **c2_, **c1_, **b1, **b2;
    double score;
};

END_JN

