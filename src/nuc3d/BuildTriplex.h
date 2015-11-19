#ifndef JIAN_NUC3D_BUILDTRIPLEX_H
#define JIAN_NUC3D_BUILDTRIPLEX_H

#include <util/std.h>
#include <util/Matr_.h>
#include <nuc2d/Seq2Tri.h>
#include <dg/DG.h>

namespace jian {
namespace nuc3d {

class BuildTriplex {
public:
    BuildTriplex();
    void operator ()(std::string sequence, nuc2d::Seq2Tri::TupleInfo ss_info, int number);
    std::vector<std::vector<int>> info_to_ct(const nuc2d::Seq2Tri::TupleInfo &info);
    void get_scaffold();

    MatrixXf _bound;
    std::vector<std::vector<int>> _ct;
    MatrixXf _scaffold;
    std::string _seq;
};

} /// namespace nuc3d
} /// namespace jian

#endif




