#ifndef JIAN_PDB_IFCHAIN_H
#define JIAN_PDB_IFCHAIN_H

#include "IFRes.h"
#include "Chain.h"

namespace jian {
namespace pdb {

class IFChain : public BasicChain<IFRes> {
public:
    IFChain() = default;

    IFChain(MolFile &pdb_file) : BasicChain<IFRes>(pdb_file) {}
 
};


} // namespace pdb
} // namespace jian

#endif

