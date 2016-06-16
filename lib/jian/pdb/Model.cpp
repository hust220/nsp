#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <map>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include "Model.hpp"
#include "MolFileParser.hpp"
#include "../utils/file.hpp"

namespace jian {

Model::Model(const std::string &f) {
    Residue residue;
    Chain chain;
    std::string file_name = file::name(f);
    std::string file_type = file::type(f);
//    std::cout << file_type << std::endl;
    MolFileParser *parser = MolFileParser::s_parsers[file_type](f);
    MolParsedLine *line = NULL, *old_line = NULL;
    name = file_name;
    while ((line = parser->line()) != NULL && (old_line == NULL || line->model_num == old_line->model_num)) {
        if (line->atom_name[0] == 'H') {
            delete line;
            continue;
        }
        if (old_line != NULL) {
            if (line->res_num != old_line->res_num || line->res_name != old_line->res_name || line->res_flag != old_line->res_flag) {
                residue.name = old_line->res_name;
                residue.num = old_line->res_num;
                chain.push_back(std::move(residue));
            }
            if (line->chain_name != old_line->chain_name) {
                chain.name = old_line->chain_name;
                this->push_back(std::move(chain));
            }
            delete old_line;
        }
        residue.push_back(Atom(line->atom_name, line->atom_num, line->x, line->y, line->z));
        old_line = line;
    }
    if (old_line != NULL) {
        residue.name = old_line->res_name;
        residue.num = old_line->res_num;
        chain.push_back(std::move(residue));
        chain.name = old_line->chain_name;
        this->push_back(std::move(chain));
        delete old_line;
        delete parser;
    }
}

std::string seq(const Model &model) {
    std::string seq;
    for (auto &&chain : model) for (auto &&res : chain) { seq += res.name; }
    return seq;
}

int num_residues(const Model &model) {
    int i = 0; for (auto &&chain : model) for (auto &&res : chain) i++;
    return i;
}

int num_atoms(const Model &model) {
    int i = 0; for (auto &&chain : model) for (auto &&res : chain) for (auto &&atom : res) i++;
    return i;
}

bool is_empty(const Model &model) {
    return num_residues(model) == 0;
}

std::ostream &operator <<(std::ostream &output, const Model &model) {
    int atom_num = 1;
    int residue_num = 1;
    output << std::fixed << std::setprecision(3);
    for (auto &&chain: model) {
        for (auto &&residue: chain) {
            for (auto &&atom: residue) {
                std::string atom_name = boost::replace_all_copy(atom.name, "*", "'");
                if (residue_num == 1 and std::set<std::string>{"P", "O1P", "O2P"}.count(atom_name)) continue;
                output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                        atom_num % atom_name % residue.name % chain.name % residue_num % 
                                        atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
                atom_num++;
            }
            residue_num++;
        }
        output << "TER" << std::endl;
    }
    return output;
}

Model RNA(const Model &model) {
    Model rna;
    static std::set<std::string> names {"A", "U", "G", "C"};
    for (auto &&chain: model) {
        Chain temp_chain; temp_chain.name = chain.name;
        for (auto &&residue: chain) {
            auto res = residue; res.name = res.name.substr(0, 1);
            if (names.count(res.name)) temp_chain.push_back(std::move(res));
        }
        if (!temp_chain.empty()) rna.push_back(temp_chain);
    }
    rna.name = model.name; rna.type = "RNA";
    return rna;
}

Model RNA(const std::string &s) {
    return RNA(Model(s));
}

Model DNA(const Model &model) {
    Model dna;
    static std::set<std::string> names {"DA", "DT", "DG", "DC"};
    for (auto &&chain: model) {
        Chain temp_chain; temp_chain.name = chain.name;
        for (auto &&residue: chain) {
            auto res = residue; res.name = res.name.substr(0, 2);
            if (names.count(res.name)) temp_chain.push_back(std::move(res));
        }
        if (!temp_chain.empty()) dna.push_back(temp_chain);
    }
    dna.name = model.name; dna.type = "DNA";
    return dna;
}

Model DNA(const std::string &s) {
    return DNA(Model(s));
}

Model R5P(const Model &model) {
    static std::map<std::string, std::set<std::string>> names {
        {"A", {"C5*", "O3*", "C1*", "N6", "C2"}},
        {"U", {"C5*", "O3*", "C1*", "O2", "O4"}},
        {"G", {"C5*", "O3*", "C1*", "O6", "N2"}},
        {"C", {"C5*", "O3*", "C1*", "O2", "N4"}}
    };
    Model m; m.type = "R5P";
    for (auto &&chain : model) {
        Chain new_chain; new_chain.name = chain.name; new_chain.type = m.type;
        for (auto &&res : chain) {
            Residue new_residue; new_residue.name = res.name;
            for (auto &&atom : res) if (names[res.name].count(atom.name)) new_residue.push_back(atom);
            if (!new_residue.empty()) new_chain.push_back(new_residue);
        }
        if (!new_chain.empty()) m.push_back(new_chain);
    }
    return m;
}

Model R5P(const std::string &s) {
    return R5P(RNA(s));
}

void write_pdb(const Model &model, const std::string &name) {
    std::ofstream ofile(name.c_str()); ofile << model; ofile.close();
}

} // namespace jian

