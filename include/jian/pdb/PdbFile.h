#ifndef JIAN_PDBFILE_H
#define JIAN_PDBFILE_H

#include "../util.h"
#include "MolFile.h"

namespace jian {

class PdbFile : public MolFile, public std::map<std::string, std::vector<std::string>> {
public:
    ///////////////////////////////////////
    void next();
    bool eof();
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

    PdbFile(std::string file_name) {
        if(file_name.size() <= 4 || file_name.substr(file_name.size() - 4, 4) != ".pdb") {
            die(std::string(file_name) + " is not an effective pdb file name");    
        }
        _name = file_name.substr(0, file_name.size() - 4);
        std::ifstream ifile(file_name.c_str());
        while (read_line(ifile));
        ifile.close();
    }

    int read_line(std::ifstream &ifile) {
        std::string line;
        while (true) {
            std::getline(ifile, line);
            if (!ifile) return 0;
            if (line.substr(0, 5) == "MODEL") {
                _model_num++;    
                continue;
            }
            if (line.substr(0, 4) != "ATOM") continue;
            if (line.substr(12, 4).find("H") != std::string::npos) continue;
            (*this)["atom_name"].push_back(boost::trim_copy(line.substr(12, 4)));
            (*this)["res_name"].push_back(boost::trim_copy(line.substr(17, 3)));
            (*this)["chain_name"].push_back(boost::trim_copy(line.substr(20, 2)));
            (*this)["atom_num"].push_back(boost::trim_copy(line.substr(6, 5)));
            (*this)["res_num"].push_back(boost::trim_copy(line.substr(22, 4)));
            (*this)["model_num"].push_back(std::to_string(_model_num));
            (*this)["x"].push_back(boost::trim_copy(line.substr(30, 8)));
            (*this)["y"].push_back(boost::trim_copy(line.substr(38, 8)));
            (*this)["z"].push_back(boost::trim_copy(line.substr(46, 8)));
            return 1;
        }
    }

};

inline void PdbFile::next() {
    if (!eof()) {
        _i++;
    }
}

inline bool PdbFile::eof() {
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

