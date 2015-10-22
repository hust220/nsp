#ifndef JIAN_NUC2D_SEQ2SS_H
#define JIAN_NUC2D_SEQ2SS_H

#include <util/util.h>

namespace jian {
namespace nuc2d {

class Seq2Ss {
public:
    std::vector<std::pair<int, int>> operator ()(std::string seq);
    int mms(int, int);
    int score(int, int);
    void backtrack(int, int);
    int get_mms(int, int);

    MatrixXf _mms;
    std::string _seq;
    std::vector<int> _types;
    int _len;
    int _min_hairpin_size = 3;
    std::vector<std::pair<int, int>> _pairs;
};

} /// namespace nuc2d
} /// namespace jian

#endif
