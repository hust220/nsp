#ifndef JIAN_GEOM_UTIL_H
#define JIAN_GEOM_UTIL_H

namespace jian {
namespace geom {

template<typename Mat1, typename Mat2> double rmsd(Mat1 &&mat1, Mat2 &&mat2) {
    SupPos sp;
    sp(mat1, mat2);
    return sp.rmsd;
}




} // namespace geom
} // namespace jian









#endif

