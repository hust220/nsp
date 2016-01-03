#ifndef JIAN_GEOM_MOVE_H
#define JIAN_GEOM_MOVE_H

#include "../util/util.h"

namespace jian {
namespace geom {

template<typename T, typename Mat> inline void move(T &&t, Mat &&mat) {
    for (int i = 0; i < 3; i++) t[i] += mat[i];
}

} // namespace geom
} // namespace jian






#endif

