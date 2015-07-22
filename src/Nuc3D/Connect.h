#ifndef CONNECT_H
#define CONNECT_H

#include "../Pdb.h"

namespace jian {

namespace nuc3d {

class Connect {
public:
//	~Connect() {
//		if (_a) delete _a;
//		if (_b) delete _b;
//	}
	Model operator ()(const Model &, const Model &, int, int);
	static void translate(Model &, const Point &);
	static void rotate(Model &, const Matrix3f &);

	Model _model_1;
	Model _model_2;
	set<string> _superposed_atoms{"C5*", "O3*", "C1*"};

}; // Connect

} /// namespace nuc3d

} /// namespace jian

#endif //CONNECT_H







