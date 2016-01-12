#ifndef JIAN_MOLFILE_H
#define JIAN_MOLFILE_H

#include "../util/util.h"

namespace jian {

class MolFile {
public:
    virtual void next() = 0;
    virtual bool eof() = 0;
    virtual double x() = 0;
    virtual double y() = 0;
    virtual double z() = 0;
    virtual std::string atom_name() = 0;
    virtual std::string atom_type() = 0;
    virtual int atom_num() = 0;
    virtual std::string res_name() = 0;
    virtual int res_num() = 0;
    virtual std::string chain_name() = 0;
    virtual int model_num() = 0;

    std::string _name;
    unsigned int _i = 0;
};

} /// namespace jian

#endif

