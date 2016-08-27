#include <iostream>
#include <iomanip>
#include <string>
#include <set>
#include <list>
#include "../utils/file.hpp"
#include "io.hpp"

namespace jian {

void mol_write_line(
    std::ostream & output,
    const int         & atom_num,
    const std::string & atom_name,
    const std::string & residue_name,
    const std::string & chain_name,
    const int         & residue_num,
    const double      & x,
    const double      & y,
    const double      & z,
    const double      & a,
    const double      & b,
    const char        & atom_name_label
) {
    output << std::fixed
           << "ATOM" 
           << std::setw(7)  << atom_num
           << "  "
           << std::left
           << std::setw(4)  << atom_name
           << std::right
           << std::setw(3)  << std::right << residue_name
           << std::setw(2)  << chain_name
           << std::setw(4)  << residue_num
           << std::setprecision(3)
           << std::setw(12) << x
           << std::setw(8)  << y
           << std::setw(8)  << z
           << std::setprecision(2)
           << std::setw(6)  << a
           << std::setw(6)  << b
           << std::setw(12) << atom_name_label
           << "  "
           << std::endl;
}

bool diff_model(const molstream &parser) {
    auto && line = parser._next_line;
    auto && old_line = parser._curr_line;
    return line->model_num != old_line->model_num;
}

bool diff_chain(const molstream &parser) {
    auto && line = parser._next_line;
    auto && old_line = parser._curr_line;
    return line->chain_name != old_line->chain_name || 
           diff_model(parser);
}

bool diff_residue(const molstream &parser) {
    auto && line = parser._next_line;
    auto && old_line = parser._curr_line;
    return line->res_num != old_line->res_num || line->res_name != old_line->res_name || line->res_flag != old_line->res_flag ||
           diff_chain(parser);
}

bool diff_atom(const molstream &parser) {
    auto && line = parser._next_line;
    auto && old_line = parser._curr_line;
    return line->atom_num != old_line->atom_num || line->atom_name != old_line->atom_name || 
           diff_residue(parser);
}

void chain_read_model(Chain &s, std::string f, std::string type) {
    molstream *parser = molstream::make(jian::file::type(f), f, type);
    do {
        (*parser) >> s;
    } while (!parser->eof() && !diff_model(*parser));
    delete parser;
}

Chain read_model_to_chain(std::string f, std::string type) {
    Chain chain;
    chain_read_model(chain, f, type);
    return chain;
}

void append_chain_to_file(const Chain &chain, const std::string &file_name, int n) {
    std::ofstream output(file_name.c_str(), std::ios::app);
    output << "MODEL " << n << std::endl;
    output << chain;
    output << "ENDMDL" << std::endl;
    output.close();
}

molstream &operator >>(molstream &parser, Atom &atom) {
    MolParsedLine *line = parser.parse_line();
    if (line != NULL) {
        atom.init(line->atom_name, line->x, line->y, line->z, line->atom_num);
    }
}

molstream &operator >>(molstream &parser, Residue &residue) {
//    std::cout << "molstream >> residue" << std::endl;
    residue.name = parser._next_line->res_name;
    residue.num = parser._next_line->res_num;
    for (int i = 0; !parser.eof(); i++) {
        if (i == 0 || !diff_residue(parser)) {
            Atom atom;
            parser >> atom;
            residue.push_back(atom);
        } else {
            break;
        }
    }
}

molstream &operator >>(molstream &parser, Chain &chain) {
//    std::cout << "molstream >> chain" << std::endl;
    chain.name = parser._next_line->chain_name;
    chain.model_name = parser.file_name;
    for (int i = 0; !parser.eof(); i++) {
        if (i == 0 || !diff_chain(parser)) {
            Residue residue;
            parser >> residue;
            if (res_is_type(residue, parser.mol_type)) {
                chain.push_back(residue);
            }
        } else {
            break;
        }
    }
}

molstream &operator >>(molstream &parser, Model &model) {
//    std::cout << "molstream >> model" << std::endl;
    model.num = parser._next_line->model_num;
    model.name = parser.file_name;
    for (int i = 0; !parser.eof(); i++) {
        if (i == 0 || !diff_model(parser)) {
            Chain chain;
            parser >> chain;
            model.push_back(chain);
        } else {
            break;
        }
    }
}

molstream &operator >>(molstream &parser, Molecule &mol) {
//    std::cout << "molstream >> mol" << std::endl;
    mol.name = parser.file_name;
    for (int i = 0; !parser.eof(); i++) {
        Model model;
        parser >> model;
        mol.push_back(model);
    }
}

std::ostream &operator <<(std::ostream &output, const Molecule &mol) {
    int atom_num = 1;
    int residue_num = 1;
    int model_num = 1;
    output << std::fixed << std::setprecision(3);
    for (auto && model : mol) {
        output << "MODEL " << model_num << std::endl;
        for (auto && chain : model) {
            for (auto &&residue: chain) {
                for (auto &&atom: residue) {
                    std::string atom_name(atom.name);
                    std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
                    mol_write_line(output, atom_num, atom_name, residue.name, chain.name, residue_num, atom[0], atom[1], atom[2], 1, 0, atom_name[0]);
                    atom_num++;
                }
                residue_num++;
            }
            output << "TER" << std::endl;
        }
        output << "ENDMDL" << std::endl;
        model_num++;
    }
}

std::ostream &operator <<(std::ostream &output, const Model &model) {
    int atom_num = 1;
    for (auto &&chain: model) {
        int residue_num = 1;
        for (auto &&residue: chain) {
            for (auto &&atom: residue) {
                std::string atom_name = atom.name;
                std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
                if (residue_num == 1 and std::set<std::string>{"P", "O1P", "O2P"}.count(atom_name)) continue;
                mol_write_line(output, atom_num, atom_name, residue.name, chain.name, residue_num, atom[0], atom[1], atom[2], 1, 0, atom_name[0]);
                atom_num++;
            }
            residue_num++;
        }
        output << "TER" << std::endl;
    }
    return output;
}

std::ostream &operator <<(std::ostream &output, const Chain &chain) {
    int atom_num = 1;
    int residue_num = 1;
    for (auto &&residue: chain) {
        for (auto &&atom: residue) {
            std::string atom_name = atom.name;
            std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
            mol_write_line(output, atom_num, atom_name, residue.name, chain.name, residue_num, atom[0], atom[1], atom[2], 1, 0, atom_name[0]);
            atom_num++;
        }
        residue_num++;
    }
    output << "TER" << std::endl;
}

std::ostream &operator <<(std::ostream &output, const Residue &residue) {
    int atom_num = 1;
    for (auto &&atom: residue) {
        std::string atom_name = atom.name;
        std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
        mol_write_line(output, atom_num, atom_name, residue.name, "X", 1, atom[0], atom[1], atom[2], 1, 0, atom_name[0]);
        atom_num++;
    }
}

std::ostream &operator <<(std::ostream &output, const Atom &atom) {
    std::string atom_name = atom.name;
    std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
    mol_write_line(output, 1, atom_name, "X", "X", 1, atom[0], atom[1], atom[2], 1, 0, atom_name[0]);
}

}

