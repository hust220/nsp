#ifndef JIAN_NUC2D_SEQ2SS_H
#define JIAN_NUC2D_SEQ2SS_H

#include <util/util.h>

namespace jian {
namespace nuc2d {

class Seq2Ss {
public:
    typedef std::pair<int, int> Pair;
    typedef std::vector<Pair> Pairs;
    typedef std::vector<int> PairList;
    typedef std::pair<PairList, double> PairInfo;
    typedef std::vector<PairInfo> InfoList;

    Seq2Ss();
    void operator ()(std::string seq);
//    double mms(int, int);
    InfoList sub_info(int, int, double);
//    double get_mms(int, int);
    InfoList best_info(int m, int n);
    PairInfo strand_info(int m, int n);
    PairInfo score(const PairInfo &info1, int m, const PairInfo &info2, int n);
    PairInfo score(const PairInfo &info, int m);
    std::pair<double, bool> pair_energy(int i, int j);
    double stack_energy(int i, int j);


    std::map<int, InfoList> _info;

    MatrixXf _mms;
    std::string _seq;
    std::vector<int> _types;
    int _len;
    int _min_hairpin_size = 4;
    double _cutoff = 0.2;
    Matrix4f _pair_energy;
    Matrix4f _stack_energy;
};

} /// namespace nuc2d
} /// namespace jian

#endif
