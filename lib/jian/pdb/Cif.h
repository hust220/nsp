#ifndef CIF_H
#define CIF_H

#include "../util/util.h"
#include "MolFile.h"

namespace jian {

class Cif : public MolFile, public std::map<std::string, std::string> {
public:
    char _ch = '\0';
    std::map<std::string, std::vector<std::string>> _loop;

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


    Cif(std::string file_name) {
        std::ifstream ifile(file_name.c_str());
        read_name(ifile);
        while (read_block(ifile));
        ifile.close();
    }

    void read_name(std::ifstream &ifile) {
        _name = read_word(ifile);
        _name = _name.substr(_name.size() - 4, 4);
        read_word(ifile);
    }

    int read_block(std::ifstream &ifile) {
        std::string key = read_word(ifile);
        if (key == "") {
            return 0;    
        } else if (key == "#") {
            return 1;
        }
        std::string value;
        if (key == "loop_") {
            read_loop(ifile);
            return 1;
        } else {
            while (true) {
                value = read_word(ifile);
                if (value == "" || value == "#") {
                    die("Cif::read_block error!");
                }
                (*this)[key] = value;
                key = read_word(ifile);
                if (key == "#") {
                    return 1;
                } else if (key == "") {
                    return 0;
                }
            }
        }
    }

    void read_loop(std::ifstream &ifile) {
        std::vector<std::string> names;
        std::string word;
        while (true) {
            word = read_word(ifile);
            if (word[0] == '_') {
                names.push_back(word);
            } else {
                while (true) {
                    for (int i = 0; i < names.size(); i++) {
                        _loop[names[i]].push_back(word);
                        word = read_word(ifile);
                        if (word == "" || word == "#") {
                            BOOST_ASSERT(i + 1 == names.size() && "Cif::read_loop error!");
                            return;
                        }
                    }
                }
            }
        }
    }

    std::string read_word(ifstream &ifile) {
        std::string word;
        A(ifile, word);
        return word;
    }

    void A(ifstream &ifile, std::string &word) {
        char c = ifile.get();
        switch (char_type(c)) {
            case 1: 
            case 2: A(ifile, word); break;
            case 3: C(ifile, word); break;
            case 4: D(ifile, word); break;
            case 5: E(ifile, word); break;
            case 6: word += c; B(ifile, word); break;
            case 7: return;
        }
    }

    void B(ifstream &ifile, std::string &word) {
        char c = ifile.get();
        switch (char_type(c)) {
            case 1: 
            case 2: return;
            case 3: C(ifile, word); break;
            case 4: D(ifile, word); break;
            case 5: E(ifile, word); break;
            case 6: word += c; B(ifile, word); break;
            case 7: return;
        }
    }

    void C(ifstream &ifile, std::string &word) {
        char c = ifile.get();
        switch (char_type(c)) {
            case 1: word += " "; C(ifile,word); break;
            case 2: C(ifile, word); break;
            case 4: 
            case 5:
            case 6: word += c; C(ifile, word); break;
            case 3: B(ifile, word); break;
            case 7: return;
        }
    }

    void D(ifstream &ifile, std::string &word) {
        char c = ifile.get();
        switch (char_type(c)) {
            case 1: word += " "; D(ifile,word); break;
            case 2: D(ifile, word); break;
            case 3: 
            case 5:
            case 6: word += c; D(ifile, word); break;
            case 4: B(ifile, word); break;
            case 7: return;
        }
    }

    void E(ifstream &ifile, std::string &word) {
        char c = ifile.get();
        switch (char_type(c)) {
            case 1: word += " "; E(ifile,word); break;
            case 3: 
            case 4: 
            case 5:
            case 6: word += c; E(ifile, word); break;
            case 2: F(ifile, word); break;
            case 7: return;
        }
    }

    void F(ifstream &ifile, std::string &word) {
        char c = ifile.get();
        switch (char_type(c)) {
            case 1: word += " "; E(ifile,word); break;
            case 3: 
            case 4: 
            case 6: word += c; E(ifile, word); break;
            case 2: F(ifile, word); break;
            case 5: return;
            case 7: return;
        }
    }

    int char_type(char c) {
        if (c == ' ' || c == '\t') { 
            return 1;
        } else if (c == '\r' || c == '\n') {
            return 2;
        } else if (c == '\'') {
            return 3;
        } else if (c == '"') {
            return 4;
        } else if (c == ';') {
            return 5;
        } else if (c == std::char_traits<char>::eof()) {
            return 7;
        } else {
            return 6;
        }
    }

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

