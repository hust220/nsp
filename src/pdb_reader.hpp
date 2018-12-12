#pragma once

#include "pdb.hpp"

namespace jian {

class PdbReader {
public:
    Residue atoms;
    Chain residues;
    Model chains;
    Molecule &models;

    PdbReader(Molecule &pdb) : models(pdb) {}

    int model_num = 0;

    struct ParsedLine {
        Str atom_name, atom_type, atom_flag, res_name, res_flag, chain_name;
        int atom_num, res_num;
        double x, y, z;
        bool is_std;
    };

    ParsedLine ol;

    void parse_line(Str line, ParsedLine &pl) {
        pl.res_name   = string_trim_c(line.substr(17, 3));
        pl.res_num    = JN_INT(string_trim_c(line.substr(22, 4)));
        pl.res_flag   = string_trim_c(line.substr(26, 1));
        pl.chain_name = string_trim_c(line.substr(20, 2));
        pl.is_std     = (!line.compare(0, 4, "ATOM"));

        pl.atom_name  = string_trim_c(line.substr(12, 4));
        pl.atom_num   = JN_INT(string_trim_c(line.substr(6, 5)));
        pl.atom_flag  = string_trim_c(line.substr(16, 1));
        pl.x          = JN_DBL(string_trim_c(line.substr(30, 8)));
        pl.y          = JN_DBL(string_trim_c(line.substr(38, 8)));
        pl.z          = JN_DBL(string_trim_c(line.substr(46, 8)));

        if (line.size() >= 78) {
            pl.atom_type  = string_trim_c(line.substr(76, 2));
        }
        else {
            pl.atom_type  = "X";
        }

//        std::replace(pl.atom_name.begin(), pl.atom_name.end(), '*', '\'');
//        if (pl.atom_name == "O1P") pl.atom_name = "OP1";
//        if (pl.atom_name == "O2P") pl.atom_name = "OP2";
        std::replace(pl.atom_name.begin(), pl.atom_name.end(), '\'', '*');
        if (pl.atom_name == "OP1") pl.atom_name = "O1P";
        if (pl.atom_name == "OP2") pl.atom_name = "O2P";
    };

    void add_residue() {
        if (!atoms.empty()) {
            Residue residue = std::move(atoms);
            residue.is_std = std::all_of(atoms.begin(), atoms.end(), [](const Atom &atom){ return atom.is_std; });
            residue.name = ol.res_name;
            residue.num = ol.res_num;
            residues.push_back(std::move(residue));
        }
    };

    void add_chain() {
        add_residue();
        if (!residues.empty()) {
            Chain chain = std::move(residues);
            chain.name = ol.chain_name;
            chains.push_back(std::move(chain));
        }
    };

    void add_atom(Str line) {
        ParsedLine pl;
        parse_line(line, pl);

        if (!atoms.empty()) {
            if (pl.res_num != ol.res_num || pl.res_name != ol.res_name || pl.res_flag != ol.res_flag || pl.chain_name != ol.chain_name) {
                add_residue();
            }
            if (pl.chain_name != ol.chain_name) {
                add_chain();
            }
        }

        if (string_starts_with(pl.atom_name, "H")) return;
        if (pl.res_name == "HOH" || pl.res_name == "H2O") return;

        // Check whether the atom has appeared before.
        if (!pl.atom_flag.empty()) {
            if (std::any_of(atoms.begin(), atoms.end(), [&pl](const Atom &atom){
                return atom.name == pl.atom_name;
            })) return;
        }

        Atom atom;
        atom[0] = pl.x;
        atom[1] = pl.y;
        atom[2] = pl.z;
        atom.name = pl.atom_name;
        atom.type = pl.atom_type;
        atom.num = pl.atom_num;
        atom.is_std = pl.is_std;

        atoms.push_back(std::move(atom));

        ol = pl;
    };

    void add_model() {
        add_chain();
        if (!chains.empty()) {
            Model model = std::move(chains);
            model.num = model_num;
            models.push_back(std::move(model));
            model_num++;
        }
    };

    void read(const std::string &fn) {
        std::ifstream ifile(fn.c_str());
        while (ifile) {
            Str line;
            std::getline(ifile, line);
            if (string_starts_with(line, "ATOM")) {
                add_atom(line);
            }
            else if (string_starts_with(line, "HETATM")) {
                add_atom(line);
            }
            else if (string_starts_with(line, "MODEL")) {
                add_model();

                auto &&v = string_tokenize(line, " ");
                if (v.size() >= 2) {
                    model_num = JN_INT(v[1])-1;
                }
            }
            else if (string_starts_with(line, "ENDMDL")) {
                add_model();
            }
            else if (string_starts_with(line, "TER")) {
                add_chain();
            }
            else if (string_trim_c(line) == "END") {
                break;
            }
            else {
                continue;
            }
        }
        ifile.close();

        add_model();

        models.name = fn;
    }
};

} // namespace jian

