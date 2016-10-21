#include <iostream>
#include <string>
#include "molstream.hpp"

namespace jian {

molstream::molstream(const std::string &f) {
	FOPEN(ifile, f);
}

molstream::~molstream() {
	FCLOSE(ifile);
    for (auto && l : gc_line) {
        delete l;
    }
}

MolParsedLine *molstream::parse_line() {
    MolParsedLine *l;
    while (true) {
        l = getline();
        if (l == NULL) {
            break;
        } else if (l->atom_name[0] == 'H') {
            delete l;
        } else {
            break;
        }
    };
    if (l != NULL) gc_line.push_back(l);
//    std::cout << "parse_line: " << l << std::endl;
    _curr_line = _next_line;
    _next_line = l;
    return _curr_line;
}

bool molstream::eof() {
    return _next_line == NULL;
}

//molstream::factory_t &molstream::parsers() {
//    static factory_t parser;
//    return parser;
//}

molstream *molstream::make(const std::string &file_type, const std::string &file_path, std::string mol_type) {
    molstream *parser = FacMolParser::create(file_type, file_path);
    parser->mol_type = mol_type;
    parser->file_name = jian::file::name(file_path);
    parser->file_type = jian::file::type(file_path);
    parser->parse_line();
    return parser;
}

//reg_molstream::reg_molstream(std::string s, molstream::creater_t f) {
//    molstream::parsers()[s] = f;
//}

} // namespace jian

