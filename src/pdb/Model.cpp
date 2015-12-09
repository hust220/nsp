#include "Model.h"

namespace jian {

Model::Model() {
}

Model::Model(MolFile &pdb_file) {
    name = pdb_file._name;
    if (!pdb_file.eof()) {
        int num = pdb_file.model_num();
        while (!pdb_file.eof() && num == pdb_file.model_num()) {
            chains.push_back(Chain(pdb_file));
        }
    }
}

Model::Model(std::string pdbfile) {
    read(pdbfile);
}

Model &Model::operator =(const Model &model) {
    name = model.name;
    chains = model.chains;
    return (*this);
}

void Model::read(string pdbfile) {
    if (pdbfile.size() > 4 && pdbfile.substr(pdbfile.size() - 4, 4) == ".pdb") {
        read_pdb(pdbfile);
    } else if (pdbfile.size() > 4 && pdbfile.substr(pdbfile.size() - 4, 4) == ".cif") {
        read_cif(pdbfile);
    } else {
        throw "JIAN::MODEL::read(std::string) error! Please give me a file ended with '.pdb' or '.cif'!";
    }
}

void Model::read_pdb(std::string pdbfile) {
    PdbFile pdb_file(pdbfile);
    (*this) = Model(pdb_file);
}

void Model::read_cif(std::string file_name) {
    Cif cif(file_name);
    (*this) = Model(cif);
}

Model::Model(vector<string> lines) {
    name = "none";
    vector<string> chain_lines;
    char chain_name = ' ';
    int atom_num = 0;

    for (auto &line: lines) {
        atom_num++;
        if (line[21] != chain_name && atom_num != 1) {
            push(Chain(chain_lines));
            chain_lines.clear();
        }
        chain_name = line[21];
        chain_lines.push_back(line);
    }
    push(Chain(chain_lines));
    chain_lines.clear();
}

int Model::empty() {
    int n = 0;
    for (auto &&chain: chains) {
        for (auto &&res: chain.residues) {
            for (auto &&atom: res.atoms) {
                n++;
            }
        }
    }
    if (n == 0) {
        return 1;
    } else {
        return 0;
    }
}

std::string Model::seq(std::string delimiter) {
    std::vector<std::string> res_names;
    for (auto &&chain: chains) {
        for (auto &&res: chain.residues) {
            res_names.push_back(res.name);
        }
    }
    return std::accumulate(res_names.begin() + 1, res_names.end(), res_names[0], [&](std::string a, std::string b){
        return a + delimiter + b;
    });
}

Chain &Model::operator [](int n) {
    return chains[n];
}

const Chain &Model::operator [](int n) const {
    return chains[n];
}

void Model::push(const Chain &chain) {
    chains.push_back(chain);
}

ostream &operator <<(ostream &output, const Model &model) {
    int atom_num = 1;
    int residue_num = 1;
    output << fixed << setprecision(3);
    for (auto &&chain: model.chains) {
        for (auto &&residue: chain.residues) {
            for (auto &&atom: residue.atoms) {
                std::string atom_name = boost::replace_all_copy(atom.name, "*", "'");
                output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % atom_num % atom_name % residue.name % chain.name % residue_num % atom.x % atom.y % atom.z % 1.00 % 0.00 % atom_name[0];
                atom_num++;
            }
            residue_num++;
        }
        output << "TER" << endl;
    }
    return output;
}

void Model::print() {
    cout << *this;
}

void Model::write(string pdbname) {
    ofstream ofile(pdbname.c_str());
    ofile << (*this);
    ofile.close();
}

vector<Chain>::iterator Model::begin() {
    return chains.begin();
}
vector<Chain>::iterator Model::end() {
    return chains.end();
}

// # Return the nth residue
// The efficient is linear to n
} /// namespace jian



