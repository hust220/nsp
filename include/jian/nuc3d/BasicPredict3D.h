#ifndef JIAN_NUC3D_BASIC_PREDICT_3D
#define JIAN_NUC3D_BASIC_PREDICT_3D

namespace jian {
namespace nuc3d {

class BasicPredict3D {
public:
    virtual Model predict() = 0;
};

} // namespace nuc3d
} // namespace jian

#endif

