#ifndef CIF_H
#define CIF_H

#include "../util.h"
#include "MolFile.h"

namespace jian {

class Cif : public MolFile, public std::map<std::string, std::string> {
public:
    char _ch = '\0';
    std::map<std::string, std::vector<std::string>> _loop;
    std::vector<std::vector<int>> _fsm = {
        // space ret ' " word ;
        { 1,  0,  2,  3,  4,  5},
        { 1,  0,  2,  3,  4,  4},
        { 2, -2, -1,  2,  2,  2},
        { 3, -2,  3, -1,  3,  3},
        {-1, -1, -2, -2,  4,  4},
        { 5,  6,  5,  5,  5,  5},
        { 5,  6,  5,  5,  5,  7},
        { 7, -1,  5,  5,  5,  5},
    };
    unsigned int _num_line = 1;

    ///////////////////////////////////////
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
        if (key == "") return 0; else if (key == "#") return 1;
        std::string value;
        if (key == "loop_") {
            read_loop(ifile);
            return 1;
        } else {
            while (true) {
                value = read_word(ifile);
                if (value == "" || value == "#") die("Cif::read_block error!");
                (*this)[key] = value;
                key = read_word(ifile);
                if (key == "#") return 1; else if (key == "") return 0;
            }
        }
    }

    void read_loop(std::ifstream &ifile) {
        std::vector<std::string> names;
        std::string word;
        while (true) {
            word = read_word(ifile);
            if (word[0] == '_') names.push_back(word);
            else {
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
        int state = 0, new_state;
        while (true) {
            char c = ifile.get();
            if (c == '\n') _num_line++;
            auto i = char_type(c);
            if (i == 6) return word;
            new_state = _fsm[state][i];
            if (new_state == -1) return strip(word);
            else if (new_state == -2) throw "jian::pdb::Cif::read_word error! Line" + 
                std::to_string(_num_line) + " in file " + _name + ".cif.";
            else if (new_state != 0 and new_state != 1) {
                word += c;    
                state = new_state;
            }
        }
    }

    std::string strip(const std::string &s) {
        if (s.size() >= 2 and (s[0] == '\'' or s[0] == '"')) {
//            std::cout << s.substr(1) << std::endl;
            return s.substr(1);
        } else if (s[0] == ';') {
//            std::cout << s.substr(1, s.size() - 2) << std::endl;
            return s.substr(1, s.size() - 2);
        } else {
//            std::cout << s << std::endl;
            return s;    
        }
    }

    int char_type(char c) {
        if (c == ' ' || c == '\t') return 0;
        else if (c == '\r' || c == '\n') return 1;
        else if (c == '\'') return 2;
        else if (c == '"') return 3;
        else if (c == ';') return 5;
        else if (c == std::char_traits<char>::eof()) return 6;
        else return 4;
    }

    bool eof() {
        return _i >= _loop["_atom_site.group_PDB"].size();
    }

    void next() {
        while (!eof()) {
            _i++;
            if (eof() || _loop["_atom_site.group_PDB"][_i] == "ATOM") break;
        }
    }

};

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

