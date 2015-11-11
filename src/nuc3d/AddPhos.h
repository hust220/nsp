#ifndef ADDPHOS_H
#define ADDPHOS_H

#include <pdb/util.h>

namespace jian {

namespace nuc3d {

class AddPhos {
public:
	void operator ()(Model &model);
};

} /// namespace nuc3d

} /// namespace jian

#endif

