#ifndef JIAN_GEOM_TRANSLATE
#define JIAN_GEOM_TRANSLATE

#include "../util/util.h"

namespace jian {
namespace geom {

template<typename T, typename Mat> 
void translate(T &&t, Mat &&mat) {
    for (int i = 0; i < 3; i++) t[i] += mat[i];
}

} // namespace geom
} // namespace jian

#endif

