#ifndef CONNECT_H
#define CONNECT_H

#include <pdb/util.h>

namespace jian {

namespace nuc3d {

class Connect {
public:
    Model operator ()(const Model &, const Model &, int, int);
    static void translate(Model &, const Point &);
    static void rotate(Model &, const Matrix3f &);

    Model _model_1;
    Model _model_2;
    int _hinge_size = 2;
    set<string> _superposed_atoms{"C5*", "O3*", "C1*"};

}; // Connect

} /// namespace nuc3d

} /// namespace jian

#endif //CONNECT_H







