#pragma once

#include "Chain.hpp"

#define EACH_INDEX_RES_HELPER(n, m, c) ({int N_RES##n = 0; EACH((CHAIN##n, N_CHAIN##n), m, EACH(RES##n, CHAIN##n, c)); N_RES##n;})
#define EACH_RES1(m, c) EACH_INDEX_RES_HELPER(1, m, c)
#define EACH_RES2(m, c) EACH_INDEX_RES_HELPER(2, m, c)
#define EACH_RES3(m, c) EACH_INDEX_RES_HELPER(3, m, c)
#define EACH_RES4(m, c) EACH_INDEX_RES_HELPER(4, m, c)
#define EACH_RES(m, c) ({int N_RES = 0; EACH((CHAIN, N_CHAIN), m, EACH(RES, CHAIN, c; N_RES++)); N_RES;})

#define EACH_ATOM_HELPER(n, m, c) ({int N_ATOM##n = 0; EACH_RES_HELPER(n, m, EACH(ATOM##n, RES##n, c)); N_ATOM##n++;})
#define EACH_ATOM1(m, c) EACH_ATOM_HELPER(1, m, c)
#define EACH_ATOM2(m, c) EACH_ATOM_HELPER(2, m, c)
#define EACH_ATOM3(m, c) EACH_ATOM_HELPER(3, m, c)
#define EACH_ATOM4(m, c) EACH_ATOM_HELPER(4, m, c)
#define EACH_ATOM(m, c) ({int N_ATOM = 0; EACH_RES(m, EACH(ATOM, RES, c)); N_ATOM++;})

namespace jian {

class Model : public std::deque<Chain> {
public:
    std::string name = "unknown";
    std::string type = "unknown";

    Model() {}

    Model(MolFile &pdb_file) {
        name = pdb_file._name;
        if (!pdb_file.eof()) {
            int num = pdb_file.model_num();
            while (!pdb_file.eof() && num == pdb_file.model_num()) {
                this->push_back(Chain(pdb_file));
            }
        }
    }

    Model(std::string pdbfile) {
        read(pdbfile);
    }

    void read(std::string pdbfile) {
        if (pdbfile.size() > 4 && pdbfile.substr(pdbfile.size() - 4, 4) == ".pdb") {
            read_pdb(pdbfile);
        } else if (pdbfile.size() > 4 && pdbfile.substr(pdbfile.size() - 4, 4) == ".cif") {
            read_cif(pdbfile);
        } else {
            throw "JIAN::MODEL::read(std::string) error! Please give me a file ended with '.pdb' or '.cif'!";
        }
    }

    void read_pdb(std::string pdbfile) {
        PdbFile pdb_file(pdbfile);
        (*this) = Model(pdb_file);
    }

    void read_cif(std::string file_name) {
        Cif cif(file_name);
        (*this) = Model(cif);
    }

    template<typename T>
    Model coarse_grained(T &&names) const {
        Model m;
        m.name = name;
        for (auto &&chain : *this) {
            m.push_back(chain.coarse_grained(names));
        }
        return m;
    }

};

template<typename T>
std::string seq(T &&model) {
    std::string seq;
    for (auto &&chain : model) for (auto &&res : chain) { seq += res.name; }
    return seq;
}

template<typename T>
int num_residues(T &&model) {
    int i = 0; for (auto &&chain : model) for (auto &&res : chain) i++;
    return i;
}

template<typename T>
int num_atoms(T &&model) {
    int i = 0; for (auto &&chain : model) for (auto &&res : chain) for (auto &&atom : res) i++;
    return i;
}

template<typename T>
inline bool is_empty(T &&model) {
    return num_residues(model) == 0;
}

template<typename T, std::enable_if_t<std::is_same<std::decay_t<T>, Model>::value, int> = 42>
std::ostream &operator <<(std::ostream &output, T &&model) {
    int atom_num = 1;
    int residue_num = 1;
    output << fixed << setprecision(3);
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
        output << "TER" << endl;
    }
    return output;
}

template<typename T, std::enable_if_t<std::is_same<Model, std::decay_t<T>>::value, int> = 42>
Model RNA(T &&model) {
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

template<typename T, std::enable_if_t<std::is_same<std::string, std::decay_t<T>>::value, int> = 42>
Model RNA(T &&s) {
    return RNA(Model(s));
}

template<typename T, std::enable_if_t<std::is_same<Model, std::decay_t<T>>::value, int> = 42>
Model DNA(T &&model) {
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

template<typename T, std::enable_if_t<std::is_same<std::string, std::decay_t<T>>::value, int> = 42>
Model DNA(T &&s) {
    return DNA(Model(s));
}

template<typename T, std::enable_if_t<std::is_same<Model, std::decay_t<T>>::value, int> = 42>
Model R5P(T &&model) {
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

template<typename T, std::enable_if_t<std::is_same<std::string, std::decay_t<T>>::value, int> = 42>
Model R5P(T &&s) {
    return R5P(RNA(s));
}

template<typename T>
void write_pdb(T &&model, const std::string &name) {
    std::ofstream ofile(name.c_str()); ofile << model; ofile.close();
}

template<typename T> 
auto sub(const Model &model, T &&t) {
    Model m;
    int res_num = 0; 
    for (auto &&chain: model) {
        Chain temp_chain; 
        temp_chain.name = chain.name;
        for (auto &&res: chain) {
            if (std::find(std::begin(t), std::end(t), res_num) != std::end(t)) temp_chain.push_back(res);
            res_num++;
        }
        if (!temp_chain.empty()) m.push_back(temp_chain);
    }
    return m;
}

template<template<typename...> class L = std::deque>
auto residues(const Model &model) {
    L<Residue> v;
    for (auto &&chain: model) for (auto &&res: chain) v.push_back(res);
    return v;
}

template<typename T, template<typename...> class L = std::deque> 
auto residues(const Model &model, T &&ls) {
    L<Residue> v;
    int i = 0; for (auto &&chain: model) for (auto &&res: chain) {
        if (std::count(std::begin(ls), std::end(ls), i)) v.push_back(res);
        i++;
    }
    return v;
}

template<typename T>
uniform_const_t<Residue, T> &residue(T &&model, int n) {
    int i = 0; for (auto &&chain: model) for (auto &&res: chain) {
        if (i == n) return res;
        i++;
    }
    throw "JIAN::MODEL::residue(int) error! Residue index out of range.";
}

template<typename T>
auto residues_to_model(T &&ls) {
    Model model; model.resize(1);
    for (auto &&res : ls) model[0].push_back(res);
    return model;
}

} // namespace jian

