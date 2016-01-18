#ifndef JIAN_PDB_MOLECULE
#define JIAN_PDB_MOLECULE

#include "Cif.h"
#include "PdbFile.h"
#include "../util/std.h"

namespace jian {
namespace pdb {

struct _Mol {}; struct _RNA {}; struct _DNA {}; struct _Pro {};
struct _AA {}; struct _PSB {};

template<typename T>
struct is_mol_file {
    enum {value = std::is_same<std::decay_t<T>, Cif>::value or std::is_same<std::decay_t<T>, PdbFile>::value};
};

class Atom : public std::array<double, 3> {
public:
    std::string _name;
    int _num;

    Atom() : std::array<double, 3>({0, 0, 0}) {}

    template<typename F, std::enable_if_t<is_mol_file<F>::value, int> = 42>
    Atom(F &&file) {
        if (!file.eof()) {
            set_name(file.atom_name());
            _num = file.atom_num();
            at(0) = file.x(); at(1) = file.y(); at(2) = file.z();
            file.next();
        }
    }

    void set_name(const std::string &name) {
        _name = name;
        std::replace(_name.begin(), _name.end(), '\'', '*');
        if (_name == "OP1") _name = "O1P"; else if (_name == "OP2") _name = "O2P";
    }

};

template<typename T, typename U>
class Residue : public std::deque<Atom> {
public:
    using mol_type = T;
    using model_type = U;

    std::string _name;
    int _num;

    template<typename F, typename V, std::enable_if_t<std::is_same<U, V>::value, int> = 42>
    Residue(const Residue<F, V> &res) : _name(res._name) {
        for (auto &&atom : res) this->push_back(atom);
    }

    template<typename F, typename V, std::enable_if_t<std::is_same<U, _PSB>::value and std::is_same<V, _AA>::value, int> = 42>
    Residue(const Residue<F, V> &res) : _name(res._name) {
        static std::set<std::string> names_base_atoms{"N1", "C2", "N3", "C4", "C5", "C6"};
        std::deque<Atom> atoms_base;
        for (auto &&atm : res) {
            Atom atom(atm);
            if (atom._name == "P") { atom._name = "P"; this->push_back(std::move(atom));    
            } else if (atom._name == "C4*") { atom._name = "S"; this->push_back(std::move(atom));
            } else if (names_base_atoms.count(atom._name)) atoms_base.push_back(std::move(atom));
        }
        this->push_back(center_atom(atoms_base));
    }

    template<typename F, typename V, std::enable_if_t<std::is_same<U, V>::value, int> = 42>
    Residue<T, U> &operator =(const Residue<F, V> &res) {
        _name = res._name;
        this->clear();
        for (auto &&atom : res) this->push_back(atom);
        return *this;
    }

    template<typename F, typename V, std::enable_if_t<std::is_same<U, _PSB>::value and std::is_same<V, _AA>::value, int> = 42>
    Residue<T, U> &operator =(const Residue<F, V> &res) {
        static std::set<std::string> names_base_atoms{"N1", "C2", "N3", "C4", "C5", "C6"};
        _name = res._name;
        this->clear();
        std::deque<Atom> atoms_base;
        for (auto &&atm : res) {
            Atom atom(atm);
            if (atom._name == "P") { atom._name = "P"; this->push_back(std::move(atom));    
            } else if (atom._name == "C4*") { atom._name = "S"; this->push_back(std::move(atom));
            } else if (names_base_atoms.count(atom._name)) atoms_base.push_back(std::move(atom));
        }
        this->push_back(center_atom(atoms_base));
        return *this;
    }

    template<typename LS>
    Atom center_atom(LS &&ls) {
        Atom atom; atom._name = "B";
        int len = ls.size();
        for (auto && atm : ls) for (int i = 0; i < 3; i++) atom[i] += atm[i];
        for (int i = 0; i < 3; i++) atom[i] /= len;
        return atom;
    }

    template<typename F, std::enable_if_t<is_mol_file<F>::value and std::is_same<U, _AA>::value, int> = 42>
    Residue(F &&file) {
        if (!file.eof()) {
            _name = format_name(file.res_name());
            _num = file.res_num();
            int model_num = file.model_num();
            std::string chain_name = file.chain_name();
            while (!file.eof() and model_num == file.model_num() and chain_name == file.chain_name() 
                   and _num == file.res_num() and _name == format_name(file.res_name())) {
                this->push_back(Atom(file));
            }
        }
    }

    template<typename F, std::enable_if_t<is_mol_file<F>::value and std::is_same<U, _PSB>::value, int> = 42>
    Residue(F &&file) : Residue(Residue<T, _AA>(file)) {}

    std::string format_name(const std::string &s) {
        std::smatch result;
        // tLeap would append a '5' after the name of first residue and '3' after the name of the last residue
        if (std::regex_match(s, result, std::regex("^(\\w+)\\d+$"))) return result[1]; else return s;
    }

    bool exists(const std::string &s) const {
        for (auto &&atom : (*this)) if (atom._name == s) return true;
        return false;
    }

    Atom &atom(const std::string &s) {
        for (auto &&atom : (*this)) if (atom._name == s) return atom;
        throw "jian::pdb::Atom::operator [](const std::string &) error!";
    }

    const Atom &atom(const std::string &s) const {
        for (auto &&atom : (*this)) if (atom._name == s) return atom;
        throw "jian::pdb::Atom::operator [](const std::string &) const error!";
    }
};

template<typename T, typename U>
class Chain : public std::deque<Residue<T, U>> {
public:
    using mol_type = T;
    using model_type = U;
    using res_type = Residue<T, U>;

    std::string _name;

    template<typename F, typename V>
    Chain(const Chain<F, V> &chain) : _name(chain._name) {
        for (auto &&res : chain) this->push_back(Residue<T, U>(res));
    }

    template<typename F, typename V>
    Chain<T, U> &operator =(const Chain<F, V> &chain) {
        _name = chain._name;
        this->clear();
        for (auto &&res : chain) this->push_back(Residue<T, U>(res));
        return *this;
    }

    template<typename F, std::enable_if_t<is_mol_file<F>::value, int> = 42>
    Chain(F &&file) {
        if (!file.eof()) {
            _name = file.chain_name();
            int model_num = file.model_num();
            while (!file.eof() and model_num == file.model_num() and _name == file.chain_name()) {
                this->push_back(Residue<T, U>(file));
            }
        }
    }
};

template<typename T, typename U>
class Model : public std::deque<Chain<T, U>> {
public:
    using mol_type = T;
    using model_type = U;
    using chain_type = Chain<T, U>;
    using res_type = Residue<T, U>;

    std::string _name;

    Model() {}

    template<typename F, typename V>
    Model(const Model<F, V> &model) : _name(model._name) {
        for (auto &&chain : model) this->push_back(Chain<T, U>(chain));
    }

    template<typename F, typename V>
    Model<T, U> &operator =(const Model<F, V> &model) {
        _name = model._name;
        this->clear();
        for (auto &&chain : model) this->push_back(Chain<T, U>(chain));
        return *this;
    }

    Model(const std::string &file_path) {
        std::string file_type = file::type(file_path);
        if (file_type == ".pdb") read_file(PdbFile(file_path));
        else if (file_type == ".cif") read_file(Cif(file_path));
        else throw "jian::pdb::Model::read(const std::string &) error! Please give me a file ended with '.pdb' or '.cif'!";
    }

    template<typename F, std::enable_if_t<is_mol_file<F>::value, int> = 42>
    Model(F &&file) {
        read_file(std::forward<T>(file));
    }

    template<typename F, std::enable_if_t<is_mol_file<F>::value, int> = 42>
    void read_file(F &&file) {
        _name = file._name;
        if (!file.eof()) {
            int num = file.model_num();
            while (!file.eof() and num == file.model_num()) this->push_back(Chain<T, U>(file));
        }
    }
};

template<typename T, typename U> 
std::ostream &operator <<(std::ostream &output, const Model<T, U> &model) {
    int atom_num = 1;
    int residue_num = 1;
    output << fixed << setprecision(3);
    for (auto &&chain: model) {
        for (auto &&residue: chain) {
            for (auto &&atom: residue) {
                std::string atom_name = boost::replace_all_copy(atom._name, "*", "'");
//                if (residue_num == 1 and std::set<std::string>{"P", "O1P", "O2P"}.count(atom_name)) continue;
                output << boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n") % 
                                        atom_num % atom_name % residue._name % chain._name % residue_num % 
                                        atom[0] % atom[1] % atom[2] % 1.00 % 0.00 % atom_name[0];
                atom_num++;
            }
            residue_num++;
        }
        output << "TER" << endl;
    }
    return output;
}

using RNA = Model<_RNA, _AA>;
using DNA = Model<_DNA, _AA>;
using Pro = Model<_Pro, _AA>;
using AA = Model<_Mol, _AA>;
using RNAPSB = Model<_RNA, _PSB>;
using DNAPSB = Model<_DNA, _PSB>;
using ProPSB = Model<_Pro, _PSB>;
using PSB = Model<_Mol, _PSB>;

} // namespace pdb
} // namespace jian

#endif

