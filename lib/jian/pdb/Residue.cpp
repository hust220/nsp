#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <algorithm>
#include "../utils/file.hpp"
#include "MolFileParser.hpp"
#include "Residue.hpp"
#include <boost/format.hpp>

namespace jian {

auto Residue::get_sort_keys() {
    static std::map<std::string, std::vector<std::string>> keys {
        {"A",{"P","O1P","O2P",
              "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
              "N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"}},
        {"U",{"P","O1P","O2P",
              "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
              "N1","C2","O2","N3","C4","O4","C5","C6"}},
        {"G",{"P","O1P","O2P",
              "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
              "N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"}},
        {"C",{"P","O1P","O2P",
              "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
              "N1","C2","O2","N3","C4","N4","C5","C6"}}
    };
    std::map<std::string, std::map<std::string, int>> sort_keys;
    int index = 0;
    for (auto && res_name : {"A", "U", "G", "C"}) {
        for (int i = 0; i < keys[res_name].size(); i++) {
            sort_keys[res_name][keys[res_name][i]] = index;
            index++;
        }
    }
    return sort_keys;
}

void Residue::sort() {
    static auto sort_keys = get_sort_keys();
    auto & keys = sort_keys;
    std::sort(this->begin(), this->end(), [&](auto &&a1, auto &&a2){
        return keys[name][a1.name] < keys[name][a2.name];
    });
}

              
std::string Residue::format_name(const std::string &s) {
    std::smatch result;
    if (std::regex_match(s, result, std::regex("^(\\w+)\\d+$"))) {
        // tLeap would append a '5' after the name of first residue
        // and '3' after the name of the last residue
        return result[1];
    } else return s;
}              
               
Atom &Residue::operator [](int n) {
    return std::deque<Atom>::operator [](n);
}              
               
const Atom &Residue::operator [](int n) const {
    return std::deque<Atom>::operator [](n);
}

Atom &Residue::operator [](const std::string &s) {
    for (auto &&atom : *this) if (atom.name == s) {return atom;}
    throw "jian::Residue::operator[] error! Not found atom!";
}

const Atom &Residue::operator [](const std::string &s) const {
    for (auto &&atom : *this) if (atom.name == s) {return atom;}
    throw "jian::Residue::operator[] error! Not found atom!";
}

Residue residue_from_file(const std::string &f) {
    Residue residue;
    std::string file_name = file::name(f);
    std::string file_type = file::type(f);
    MolFileParser *parser = MolFileParser::s_parsers[file_type](f);
    MolParsedLine *line = NULL, *old_line = NULL;
    while ((line = parser->line()) != NULL) {
        if (old_line != NULL) {
            if (line->res_num != old_line->res_num || line->res_name != old_line->res_name || line->res_flag != old_line->res_flag) {
                break;
            }
            delete old_line;
        }
        residue.push_back(Atom(line->atom_name, line->atom_num, line->x, line->y, line->z));
        old_line = line;
    }
    if (old_line != NULL) {
        residue.name = old_line->res_name;
        residue.num = old_line->res_num;
        delete old_line;
        delete parser;
    }
    return residue;
}

void residue_to_file(const Residue &residue, const std::string &file_name) {
    std::ofstream output(file_name.c_str());
    int atom_num = 1;
    output << std::fixed << std::setprecision(3);
    for (auto &&atom: residue) {
        std::string atom_name(atom.name);
        std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
        output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                atom_num % atom_name % residue.name % "X" % 1 % 
                                atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
        atom_num++;
    }
    output << "TER" << std::endl;
    output.close();
}

} // namespace jian

