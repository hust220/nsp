#ifndef JIAN_NUC3D_PREDICT_H
#define JIAN_NUC3D_PREDICT_H

#include <util/std.h>

namespace jian {
namespace nuc3d {

class Predict : public virtual JobInf, public virtual LM2, public virtual Assemble {
public:
    Predict(Par);
};

} // namespace nuc3d
} // namespace jian

#endif




