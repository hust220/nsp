#ifndef CIF_H
#define CIF_H

#include <util/util.h>
#include "MolFile.h"

namespace jian {

class Cif : public MolFile, public std::map<std::string, std::string> {
public:
    Cif(std::string file_name);
    void read_name(std::ifstream &ifile);
    int read_block(std::ifstream &ifile);
    void read_loop(std::ifstream &ifile);

    ///////////////////////////////////////
    /// read_word
    std::string read_word(ifstream &ifile);
    void A(ifstream &ifile, std::string &word);
    void B(ifstream &ifile, std::string &word);
    void C(ifstream &ifile, std::string &word);
    void D(ifstream &ifile, std::string &word);
    void E(ifstream &ifile, std::string &word);
    void F(ifstream &ifile, std::string &word);
    int char_type(char c);
    char _ch = '\0';
    ///////////////////////////////////////

    ///////////////////////////////////////
    void next();
    int eof();
    double x();
    double y();
    double z();
    std::string atom_name();
    std::string atom_type();
    int atom_num();
    std::string res_name();
    int res_num();
    std::string chain_name();
    int chain_num();
    int model_num();
    ///////////////////////////////////////

    std::map<std::string, std::vector<std::string>> _loop;
};

inline void Cif::next() {
    while (!eof()) {
        _i++;
        if (eof() || _loop["_atom_site.group_PDB"][_i] == "ATOM")
            break;
    }
}

inline int Cif::eof() {
    if (_i >= _loop["_atom_site.group_PDB"].size()) {
        return 1;
    } else {
        return 0;
    }
}

inline double Cif::x() {
    return std::stod(_loop["_atom_site.Cartn_x"][_i]);
}

inline double Cif::y() {
    return std::stod(_loop["_atom_site.Cartn_y"][_i]);
}

inline double Cif::z() {
    return std::stod(_loop["_atom_site.Cartn_z"][_i]);
}

inline std::string Cif::atom_name() {
    return _loop["_atom_site.label_atom_id"][_i];
}

inline std::string Cif::atom_type() {
    return _loop["_atom_site.type_symbol"][_i];
}

inline int Cif::atom_num() {
    return std::stoi(_loop["_atom_site.id"][_i]);
}

inline std::string Cif::res_name() {
    return _loop["_atom_site.label_comp_id"][_i];
}

inline int Cif::res_num() {
    std::string str = _loop["_atom_site.label_seq_id"][_i];
    if (str == ".") 
        return -1;
    else 
        return std::stoi(str);
}

inline std::string Cif::chain_name() {
    return _loop["_atom_site.label_asym_id"][_i];
}

inline int Cif::chain_num() {
    return std::stoi(_loop["_atom_site.label_entity_id"][_i]);
}

inline int Cif::model_num() {
    return std::stoi(_loop["_atom_site.pdbx_PDB_model_num"][_i]);
}

} /// namespace jian

#endif

