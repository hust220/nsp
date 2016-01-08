#ifndef JIAN_NUC3D_SCAFFOLDTOALLATOM
#define JIAN_NUC3D_SCAFFOLDTOALLATOM

namespace jian {
namespace nuc3d {

template<typename ModelType>
class ScaffoldToAllAtom {
public:
    template<typename MatType>
    ModelType operator ()(MatType &&mat) {
        return scaffold_to_all_atom(mat);
    }

    template<typename MatType>
    ModelType scaffold_to_all_atom(MatType &&mat) {
        return ModelType();
    }
};





} // namespace nuc3d
} // namespace jian

#endif

