#ifndef JIAN_PDBFILE_H
#define JIAN_PDBFILE_H

#include "../Utils.h"
#include "MolFile.h"

namespace jian {

class PdbFile : public MolFile, public std::map<std::string, std::vector<std::string>> {
public:
    PdbFile(std::string file_name);
    int read_line(std::ifstream &ifile);

    ///////////////////////////////////////
    void next();
    int eof();
    double x();
    double y();
    double z();
    std::string atom_name();
    std::string atom_type();
    std::string res_name();
    std::string chain_name();
    int atom_num();
    int res_num();
    int model_num();
    ///////////////////////////////////////

    int _model_num = 0;
};

inline void PdbFile::next() {
    if (!eof()) {
        _i++;
    }
}

inline int PdbFile::eof() {
    if (_i >= (*this)["atom_num"].size()) {
        return 1;
    } else {
        return 0;
    }
}

inline double PdbFile::x() {
    return std::stod((*this)["x"][_i]);
}

inline double PdbFile::y() {
    return std::stod((*this)["y"][_i]);
}

inline double PdbFile::z() {
    return std::stod((*this)["z"][_i]);
}

inline std::string PdbFile::atom_name() {
    return (*this)["atom_name"][_i];
}

inline std::string PdbFile::atom_type() {
    return (*this)["atom_name"][_i].substr(0, 1);
}

inline int PdbFile::atom_num() {
    return std::stoi((*this)["atom_num"][_i]);
}

inline std::string PdbFile::res_name() {
    return (*this)["res_name"][_i];
}

inline int PdbFile::res_num() {
    std::string str = (*this)["res_num"][_i];
    return std::stoi(str);
}

inline std::string PdbFile::chain_name() {
    return (*this)["chain_name"][_i];
}

inline int PdbFile::model_num() {
    return std::stoi((*this)["model_num"][_i]);
}

} /// namespace jian

#endif

