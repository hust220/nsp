#ifndef JIAN_PDB_IFMODEL_H
#define JIAN_PDB_IFMODEL_H

#include "IFChain.h"
#include "Model.h"

namespace jian {
namespace pdb {

class IFModel : public BasicModel<IFChain> {
public:
    IFModel() = default;

    IFModel(const std::string &pdb_file) : BasicModel<IFChain>(pdb_file) {}
};


} // namespace pdb
} // namespace jian

#endif

