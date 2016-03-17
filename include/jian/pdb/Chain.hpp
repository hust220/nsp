#pragma once

#include "Residue.hpp"

namespace jian {

class Chain : public std::deque<Residue> {
public:
    std::string name = "X";
    std::string type = "unknown";

    Chain() {}

    Chain(MolFile &pdb_file) {
        if (!pdb_file.eof()) {
            name = pdb_file.chain_name();
            int model_num = pdb_file.model_num();
            while (!pdb_file.eof() && model_num == pdb_file.model_num() && name == pdb_file.chain_name()) {
                this->push_back(Residue(pdb_file));
            }
        }
    }

};

template<typename T>
Chain residues_from_file(T &&file_name) {
    Chain residues;
    if (file::type(file_name) == ".pdb") {
        PdbFile pdb_file(file_name);
        if (!pdb_file.eof()) {
            while (!pdb_file.eof()) {
                residues.push_back(Residue(pdb_file));
            }
        }
    } else if (file_name.size() > 4 && file_name.substr(file_name.size() - 4, 4) == ".cif") {
        Cif cif(file_name);
        if (!cif.eof()) {
            while (!cif.eof()) {
                residues.push_back(Residue(cif));
            }
        }
    } else {
        throw "JIAN::RESIDUE::get_residues_from_file(std::string) error! Please give me a file ended with '.pdb' or '.cif'!";
    }
    return residues;
}

template<typename T, typename U>
void residues_to_file(T &&chain, U &&file_name) {
SEELN("Reidues To File ", file_name);
    std::ofstream output(file_name.c_str());
    int atom_num = 1;
    int residue_num = 1;
    output << fixed << setprecision(3);
    for (auto &&residue: chain) {
        for (auto &&atom: residue) {
            std::string atom_name = boost::replace_all_copy(atom.name, "*", "'");
            if (residue_num == 1 and std::set<std::string>{"P", "O1P", "O2P"}.count(atom_name)) continue;
            output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                    atom_num % atom_name % residue.name % "X" % residue_num % 
                                    atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
            atom_num++;
        }
        residue_num++;
    }
    output << "TER" << endl;
    output.close();
SEELN("Reidues To File Done.");
}

} // namespace jian

