#ifndef JIAN_NUC3D_JOBINF_H
#define JIAN_NUC3D_JOBINF_H

#include <util/std.h>
#include <util/Par.h>
#include <nuc2d/util.h>

namespace jian {
namespace nuc3d {

class JobInf {
public:
    JobInf(Par);

    std::string _name = "assemble";
    std::string _seq;
    std::string _ss;
    std::string _lib;
    std::string _family = "other";
    std::string _type = "RNA";
    std::string _constraints;
    int _hinge = 2;
    int _num = 1;
    int _is_test = 0; // Is this a test case or not?
    std::string _method = "assemble";
    std::string _native;
};

} // namespace nuc3d
} // namespace jian

#endif




