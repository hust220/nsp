#include "PdbFile.h"

namespace jian {
   
PdbFile::PdbFile(std::string file_name) {
    if(file_name.size() <= 4 || file_name.substr(file_name.size() - 4, 4) != ".pdb") {
        die(std::string(file_name) + " is not an effective pdb file name");    
    }
    _name = file_name.substr(0, file_name.size() - 4);
    std::ifstream ifile(file_name.c_str());
    while (read_line(ifile));
    ifile.close();
}

int PdbFile::read_line(std::ifstream &ifile) {
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

} /// namespace jian

