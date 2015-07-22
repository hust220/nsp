#ifndef JIAN_PDB_MANIPULATE_H
#define JIAN_PDB_MANIPULATE_H

#include "Pdb.h"
#include "../Geom.h"

namespace jian {
    
namespace pdb {
    
Vector3f normal_vector(const Residue &);
Vector3f pivot(const Residue &res);

} /// namespace pdb

} /// namespace jian



#endif

