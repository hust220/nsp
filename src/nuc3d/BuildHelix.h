#ifndef JIAN_NUC3D_BUILDHELIX
#define JIAN_NUC3D_BUILDHELIX

#include <nuc2d/util.h>
#include "Connect.h"

namespace jian {
namespace nuc3d {

class BuildHelix {
public:
    BuildHelix ();
    Model operator ()(nuc2d::helix *h);
    Model operator ()(std::string seq);

    std::string _lib;
    std::string _type = "RNA";
};

} /// namespace nuc3d
} /// namespace jian

#endif

