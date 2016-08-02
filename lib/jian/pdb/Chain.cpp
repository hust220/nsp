#include <mutex>
#include <iomanip>
#include <set>
#include <boost/format.hpp>
#include "Chain.hpp"
#include "PdbFileParser.hpp"
#include "../utils/Debug.hpp"
#include "../utils/file.hpp"

namespace jian {

static std::mutex mt;

std::ostream &operator <<(std::ostream &output, const Chain &chain) {
    int atom_num = 1;
    int residue_num = 1;
    output << std::fixed << std::setprecision(3);
    for (auto &&residue: chain) {
        for (auto &&atom: residue) {
            std::string atom_name(atom.name);
            std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
            output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                    atom_num % atom_name % residue.name % chain.name % residue_num % 
                                    atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
            atom_num++;
        }
        residue_num++;
    }
    output << "TER" << std::endl;
}

Chain residues_from_file(const std::string &f) {
    std::lock_guard<std::mutex> gd(mt);
    Residue residue;
    Chain chain;
    std::string file_name = file::name(f);
    std::string file_type = file::type(f);
//std::cout << "pasers: " << MolFileParser::s_parsers.size() <<std::endl;
//    MolFileParser *parser = MolFileParser::s_parsers[file_type](f);
    MolFileParser *parser = new PdbFileParser(f);
    MolParsedLine *line = NULL, *old_line = NULL;
    chain.name = "X";
    chain.model_name = file_name;
    while ((line = parser->line()) != NULL) {
        if (old_line != NULL) {
            if (line->res_num != old_line->res_num || line->res_name != old_line->res_name || line->res_flag != old_line->res_flag) {
                residue.name = old_line->res_name;
                residue.num = old_line->res_num;
                chain.push_back(std::move(residue));
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
        delete old_line;
        delete parser;
    }
    return chain;
}

void residues_to_file(const Chain &chain, const std::string &file_name) {
//    Debug::println("Reidues To File ", file_name);
    std::ofstream output(file_name.c_str());
    int atom_num = 1;
    int residue_num = 1;
    output << std::fixed << std::setprecision(3);
    for (auto &&residue: chain) {
        for (auto &&atom: residue) {
            std::string atom_name(atom.name);
            std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
//            std::string atom_name = boost::replace_all_copy(atom.name, "*", "'");
//            if (residue_num == 1 and std::set<std::string>{"P", "O1P", "O2P"}.count(atom_name)) continue;
            output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                    atom_num % atom_name % residue.name % "X" % residue_num % 
                                    atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
            atom_num++;
        }
        residue_num++;
    }
    output << "TER" << std::endl;
    output.close();
//    Debug::println("Reidues To File Done.");
}

void append_chain_to_file(const Chain &chain, const std::string &file_name, int n) {
    std::ofstream output(file_name.c_str(), std::ios::app);
    int atom_num = 1;
    int residue_num = 1;
    output << "MODEL " << n << std::endl;
    output << std::fixed << std::setprecision(3);
    for (auto &&residue: chain) {
        for (auto &&atom: residue) {
            std::string atom_name(atom.name);
            std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
//            std::string atom_name = boost::replace_all_copy(atom.name, "*", "'");
            if (residue_num == 1 and std::set<std::string>{"P", "O1P", "O2P"}.count(atom_name)) continue;
            output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                    atom_num % atom_name % residue.name % "X" % residue_num % 
                                    atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
            atom_num++;
        }
        residue_num++;
    }
    output << "TER" << std::endl;
    output << "ENDMDL" << std::endl;
    output.close();
}

int num_atoms(const Chain &chain) {
    int n = 0;
    for (auto && res : chain) {
        for (auto && atom : res) {
            n++;
        }
    }
    return n;
}

} // namespace jian

