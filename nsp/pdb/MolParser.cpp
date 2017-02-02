#include <iostream>
#include <string>
#include "MolParser.hpp"

BEGIN_JN

MolParser::MolParser(const S &f) {
	FOPEN(ifile, f);
}

MolParser::~MolParser() {
	FCLOSE(ifile);
    for (auto && l : gc_line) {
        delete l;
    }
}

MolParsedLine *MolParser::parse_line() {
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
	if (l != NULL && l->model_num == 0) {
		if (_next_line == NULL) {
			l->model_num = 1;
		}
		else {
			l->model_num = _next_line->model_num;
		}
	}
    _curr_line = _next_line;
    _next_line = l;
    return _curr_line;
}

bool MolParser::eof() {
    return _next_line == NULL;
}

MolParser *MolParser::make(const S &file_type, const S &file_path, S mol_type) {
    MolParser *parser = FacMolParser::create(file_type, file_path);
    parser->mol_type = mol_type;
    parser->file_name = jian::file::name(file_path);
    parser->file_type = jian::file::type(file_path);
    parser->parse_line();
    return parser;
}

END_JN

