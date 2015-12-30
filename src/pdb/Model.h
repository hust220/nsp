#ifndef Model_H
#define Model_H

#include "Chain.h"

namespace jian {

template<typename ChainType>
class BasicModel {
public:
    std::string name = "none";
    std::vector<ChainType> chains;

    BasicModel() {}

    BasicModel(BasicModel<ChainType> *model) : name(model->name), chains(model->chains) {}

    BasicModel(MolFile &pdb_file) {
        name = pdb_file._name;
        if (!pdb_file.eof()) {
            int num = pdb_file.model_num();
            while (!pdb_file.eof() && num == pdb_file.model_num()) {
                chains.push_back(ChainType(pdb_file));
            }
        }
    }

    BasicModel(std::string pdbfile) {
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
        (*this) = BasicModel<ChainType>(pdb_file);
    }

    void read_cif(std::string file_name) {
        Cif cif(file_name);
        (*this) = BasicModel<ChainType>(cif);
    }

    BasicModel(std::vector<std::string> lines) {
        name = "none";
        std::vector<std::string> chain_lines;
        char chain_name = ' ';
        int atom_num = 0;

        for (auto &&line: lines) {
            atom_num++;
            if (line[21] != chain_name && atom_num != 1) {
                push(ChainType(chain_lines));
                chain_lines.clear();
            }
            chain_name = line[21];
            chain_lines.push_back(line);
        }
        push(ChainType(chain_lines));
        chain_lines.clear();
    }

    template<typename Fn> void each_res(Fn &&f) {
        int res_num = 0;
        for (auto && chain : chains) for (auto && residue : chain.residues) {
            f(residue, res_num);    
            res_num++;
        }
    }

    template<typename Fn> void each_res(Fn &&f) const {
        int res_num = 0;
        for (auto && chain : chains) for (auto && residue : chain.residues) {
            f(residue, res_num);    
            res_num++;
        }
    }

    bool empty() const {
        return res_nums() == 0;
    }

    virtual std::string seq(std::string delimiter = "") const {
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

    ChainType &operator [](int n) {
        return chains[n];
    }

    const ChainType &operator [](int n) const {
        return chains[n];
    }

    Residue &residue(int n) {
        int res_num = 0;
        for (auto &&chain: chains) for (auto &&res: chain.residues) {
                if (res_num == n) return res;
                res_num++;
        }
        throw "JIAN::MODEL::residue(int) error! Residue index out of range.";
    }

    const Residue &residue(int n) const {
        int res_num = 0;
        for (auto &&chain: chains) for (auto &&res: chain.residues) {
                if (res_num == n) return res;
                res_num++;
        }
        throw "JIAN::MODEL::residue(int) error! Residue index out of range.";
    }

    void push(const ChainType &chain) {
        chains.push_back(chain);
    }

    void print() const {
        cout << *this;
    }

    void write(string pdbname) const {
        ofstream ofile(pdbname.c_str());
        ofile << (*this);
        ofile.close();
    }

    typename std::vector<ChainType>::iterator begin() {
        return chains.begin();
    }

    typename std::vector<ChainType>::iterator end() {
        return chains.end();
    }

    int res_nums() const {
        int res_num = 0;
        for (auto &chain: chains) for (auto &residue: chain.residues) res_num++;
        return res_num;
    }

    int atom_nums() const {
        int atom_num = 0;
        for (auto &chain: chains) for (auto &residue: chain.residues) for (auto &atom: residue.atoms) atom_num++;
        return atom_num;
    }


    template<typename T> BasicModel<ChainType> sub(T &&t) const {
        BasicModel<ChainType> model;
        int res_num = 0;
        for (auto &&chain: chains) {
            ChainType temp_chain;
            temp_chain.name = chain.name;
            for (auto &&res: chain.residues) {
                if (std::count(std::begin(t), std::end(t), res_num)) {
                    temp_chain.residues.push_back(res);
                }
                res_num++;
            }
            if (!temp_chain.residues.empty())
                model.chains.push_back(temp_chain);
        }
        return model;
    }

    template<typename List> std::deque<Residue> residues(List &&list) const {
        std::deque<Residue> vec;
        int res_num = 0;
        for (auto &&chain: chains) {
            for (auto &&res: chain.residues) {
                if (std::count(std::begin(list), std::end(list), res_num)) vec.push_back(res);
                res_num++;
            }
        }
        return vec;
    }

    std::deque<Residue> residues() const {
        std::deque<Residue> vec;
        for (auto &&chain: chains) for (auto &&res: chain.residues) vec.push_back(res);
        return vec;
    }


};

typedef BasicModel<Chain> Model;

template<typename ChainType> 
std::ostream &operator <<(std::ostream &output, const BasicModel<ChainType> &model) {
    int atom_num = 1;
    int residue_num = 1;
    output << fixed << setprecision(3);
    for (auto &&chain: model.chains) {
        for (auto &&residue: chain.residues) {
            for (auto &&atom: residue.atoms) {
                std::string atom_name = boost::replace_all_copy(atom.name, "*", "'");
                output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                        atom_num % atom_name % residue.name % chain.name % residue_num % 
                                        atom.x % atom.y % atom.z % 1.00 % 0.00 % atom_name[0];
                atom_num++;
            }
            residue_num++;
        }
        output << "TER" << endl;
    }
    return output;
}

} /// namespace jian

#endif // Model_H

