#include <algorithm>
#include <iostream>
#include <string>
#include "PdbFileParser.hpp"
#include "../utils/string.hpp"

namespace jian {

reg_molstream reg_pdb_file_parser("pdb", [](const std::string &f)->molstream * {
    return new PdbFileParser(f);
});

PdbFileParser::PdbFileParser(const std::string &f) : molstream(f) {
}

MolParsedLine *PdbFileParser::getline() {
    thread_local static std::vector<std::string> v {"A5", "U5", "G5", "C5", "T5", "A3", "U3", "G3", "C3", "T3"};
    std::string line;
    std::vector<std::string> arr;
    int model_num = 1;
    while (std::getline(ifile, line)) {
        if (line.compare(0, 4, "ATOM") == 0 ||
            line.compare(0, 6, "HETATM") == 0) {
            MolParsedLine *parsed_line = new MolParsedLine;
            parsed_line->model_num = model_num;
            parsed_line->atom_name = trim_copy(line.substr(12, 4));
            parsed_line->res_name = trim_copy(line.substr(17, 3));
            if (std::find(v.begin(), v.end(), parsed_line->res_name) != v.end()) {
                parsed_line->res_name = parsed_line->res_name.substr(0, 1);
            }
            parsed_line->chain_name = trim_copy(line.substr(20, 2));
            parsed_line->atom_num = std::stoi(trim_copy(line.substr(6, 5)));
            parsed_line->res_num = std::stoi(trim_copy(line.substr(22, 4)));
            parsed_line->res_flag = trim_copy(line.substr(26, 1));
            parsed_line->x = std::stod(trim_copy(line.substr(30, 8)));
            parsed_line->y = std::stod(trim_copy(line.substr(38, 8)));
            parsed_line->z = std::stod(trim_copy(line.substr(46, 8)));
            return parsed_line;
        } else if (line.compare(0, 5, "MODEL") == 0) {
            tokenize(line, arr, " ");
            model_num = std::stoi(arr[1]);
        }
    }
    return NULL;
}

} // namespace jian

