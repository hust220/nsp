#pragma once

#include <array>
#include <string>

namespace jian {

class Chain;

class DistAnal {
public:
    using Point = std::array<double, 4>;

    double score = 0;
    double interval;
    int cutoff;
    int bins;
    double penalty = 0;

    DistAnal(double = 0.5, int = 20);
    ~DistAnal();

    void read_mol(const Chain &);
    void train();
    DistAnal &run(const Chain &);

    void read_parm(std::string);
    int res_type(std::string name);
    void free_freq(double *f);

private:
    int *num = NULL;
    int *type = NULL;
    int *ntLen = NULL;
    Point **list = NULL;
    int len = 0;
    double *freq = NULL;

};

}

