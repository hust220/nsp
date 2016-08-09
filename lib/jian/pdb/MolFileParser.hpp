#pragma once

#include <map>
#include <string>
#include <functional>
#include <fstream>

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


class MolFileParser {
public:
    MolFileParser(const std::string &f) {
        ifile.open(f);
        if (!ifile) {
            throw std::string("Open file '") + f + "' failed!";
        }
    }

    ~MolFileParser() {
        ifile.close();
    }

    virtual MolParsedLine *line() = 0;

    using parser_creater_t = std::function<MolFileParser *(const std::string &)>;
    using factory_t = std::map<std::string, parser_creater_t>;

    static factory_t &parsers() {
        static factory_t parser;
        return parser;
    }

    static MolFileParser *make(const std::string &type, const std::string &par) {
        return parsers().at(type)(par);
    }


protected:
    std::ifstream ifile;
};

class RegisterMolFileParser {
public:
    template<typename F>
    RegisterMolFileParser(const std::string &s, F &&f) {
        MolFileParser::parsers()[s] = f;
    }
};

} // namespace jian

