#ifndef JIAN_NUC2D_SEQ2SS_H
#define JIAN_NUC2D_SEQ2SS_H

#include <util/util.h>

namespace jian {
namespace nuc2d {

class Seq2Ss {
public:
    typedef std::pair<int, int> Pair;
    typedef std::vector<Pair> Pairs;
    typedef std::vector<Pair> PairList;
    typedef std::vector<std::pair<Pairs, double>> PairLists;

    void operator ()(std::string seq);
    double mms(int, int);
    double score(int, int);
    PairLists backtrack(int, int, double);
    double get_mms(int, int);

    MatrixXf _mms;
    std::string _seq;
    std::vector<int> _types;
    int _len;
    int _min_hairpin_size = 4;
    int _cutoff = 1;
};

} /// namespace nuc2d
} /// namespace jian

#endif
