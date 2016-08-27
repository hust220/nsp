#pragma once

#include <map>
#include <string>
#include <functional>
#include <fstream>
#include <list>
#include "../utils/file.hpp"

namespace jian {

struct MolParsedLine {
    std::string atom_type;
    std::string atom_name;
    std::string res_name;
    std::string chain_name;
    std::string res_flag;
    int atom_num;
    int res_num;
    int chain_num;
    int model_num;
    double x, y, z;
};


class molstream {
public:
    using creater_t = std::function<molstream *(const std::string &)>;
    using factory_t = std::map<std::string, creater_t>;

    MolParsedLine *_curr_line = NULL;
    MolParsedLine *_next_line = NULL;
    std::list<MolParsedLine *> gc_line;
    std::string file_name;
    std::string file_type;
    std::string mol_type = "";
    std::ifstream ifile;

    molstream(const std::string &f);

    ~molstream();

    MolParsedLine *parse_line();

    bool eof();

    virtual MolParsedLine *getline() = 0;

    static factory_t &parsers();

    static molstream *make(const std::string &file_type, const std::string &file_path, std::string mol_type);

};

class reg_molstream {
public:
    reg_molstream(std::string s, molstream::creater_t f);
};

} // namespace jian

