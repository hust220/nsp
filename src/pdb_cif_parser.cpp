#include "log.hpp"
#include "pdb_cif_parser.hpp"

namespace jian {

REG_MOL_PARSER("cif", CifFileParser)

class Cif : public std::map<std::string, std::string> {
public:
    S _name;
    int _i {0};
    char _ch {'\0'};
    std::map<std::string, std::vector<std::string>> _loop;
    std::vector<std::vector<int>> _fsm {
        // space return ' " word ;
        { 1,  0,  2,  3,  4,  5},
        { 1,  0,  2,  3,  4,  4},
        { 2, -2, -1,  2,  2,  2},
        { 3, -2,  3, -1,  3,  3},
        {-1, -1, -2, -2,  4,  4},
        { 5,  6,  5,  5,  5,  5},
        { 5,  6,  5,  5,  5,  7},
        { 7, -1,  5,  5,  5,  5}
    };
    int _num_line {1};

    Cif(std::ifstream &ifile) {
        read_name(ifile);
        while (read_block(ifile));
    }

    void read_name(std::ifstream &ifile) {
        _name = read_word(ifile);
//LOGI << "name: " << _name << std::endl;
        _name = _name.substr(_name.size() - 4, 4);
        read_word(ifile);
    }

    int read_block(std::ifstream &ifile) {
        S key = read_word(ifile);
        if (key == "") return 0; else if (key == "#") return 1;
        S value;
        if (key == "loop_") {
            read_loop(ifile);
            return 1;
        } else {
            while (true) {
                value = read_word(ifile);
                if (value == "" || value == "#") throw "read_block error!";
                (*this)[key] = value;
                key = read_word(ifile);
                if (key == "#") return 1; else if (key == "") return 0;
            }
        }
    }

    void read_loop(std::ifstream &ifile) {
        std::vector<std::string> names;
        S word;
        while (true) {
            word = read_word(ifile);
            if (word[0] == '_') names.push_back(word);
            else {
                while (true) {
                    for (int i = 0; i < names.size(); i++) {
                        _loop[names[i]].push_back(word);
                        word = read_word(ifile);
                        if (word == "" || word == "#") {
                            if (i + 1 != names.size()) {
                                throw "read_loop error";
                            }
                            return;
                        }
                    }
                }
            }
        }
    }

    S read_word(std::ifstream &ifile) {
        S word;
        int state = 0, new_state;
        while (true) {
            char c = ifile.get();
//            LOGI << "new line: " << c << ' ';
            if (c == '\n') _num_line++;
            auto i = char_type(c);
//            LOGI << i << ' ';
            if (i == 6) break;;
            new_state = _fsm[state][i];
//            LOGI << new_state << std::endl;
            if (new_state == -1) break;
            else if (new_state == -2) throw "jian::pdb::read_word error! Line" + 
                std::to_string(_num_line) + " in file " + _name + ".cif.";
            else if (new_state != 0 && new_state != 1) {
                word += c;    
                state = new_state;
            }
        }
//        LOGI << word << std::endl;
        return strip(word);
    }

    S strip(const S &s) {
        if (s.size() >= 2 && (s[0] == '\'' || s[0] == '"')) {
            return s.substr(1);
        } else if (s[0] == ';') {
            return s.substr(1, s.size() - 2);
        } else {
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

    MolParsedLine *line() {
        if (!eof()) {
            MolParsedLine *line = new MolParsedLine;
//std::cout << "MolParsedLine: " << line << std::endl;
            line->atom_type = _loop["_atom_site.type_symbol"][_i];
            line->atom_name = _loop["_atom_site.label_atom_id"][_i];
            line->res_name = _loop["_atom_site.label_comp_id"][_i];
            line->chain_name = _loop["_atom_site.label_asym_id"][_i];
            line->res_flag = "";
            line->atom_num = std::stoi(_loop["_atom_site.id"][_i]);
            line->res_num = res_num();
            line->chain_num = std::stoi(_loop["_atom_site.label_entity_id"][_i]);
            line->model_num = std::stoi(_loop["_atom_site.pdbx_PDB_model_num"][_i]);
            line->x = std::stod(_loop["_atom_site.Cartn_x"][_i]);
            line->y = std::stod(_loop["_atom_site.Cartn_y"][_i]);
            line->z = std::stod(_loop["_atom_site.Cartn_z"][_i]);
            next();
            return line;
        } else {
            return NULL;
        }
    }

    int res_num() {
        S str = _loop["_atom_site.label_seq_id"][_i];
        if (str == ".") 
            return -1;
        else 
            return std::stoi(str);
    }

};

CifFileParser::CifFileParser(const S &f) : MolParser(f) {
    _cif = new Cif(ifile);
}

CifFileParser::~CifFileParser() {
    delete _cif;
}

MolParsedLine *CifFileParser::getline() {
    return _cif->line();
}

}

