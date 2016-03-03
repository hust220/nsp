#ifndef JIAN_PDB_PDB
#define JIAN_PDB_PDB

#include "Model.h"

namespace jian {

template<typename ModelType>
class BasicPdb {
public:
    std::string name = "none";
    std::vector<ModelType> models;

    BasicPdb() {}

    typename std::vector<ModelType>::iterator begin() {
        return models.begin();
    }
    typename std::vector<ModelType>::iterator end() {
        return models.end();
    }

    BasicPdb(MolFile &pdb_file) {
        name = pdb_file._name;   
        while (!pdb_file.eof()) {
            models.push_back(ModelType(pdb_file));
        }
    }

    BasicPdb(std::string file_name) {
        read(file_name);
    }

    void read(std::string file_name) {
        if (file_name.size() > 4 && file_name.substr(file_name.size() - 4, 4) == ".pdb") {
            PdbFile pdb(file_name);
            (*this) = BasicPdb<ModelType>(pdb);
        } else if (file_name.size() > 4 && file_name.substr(file_name.size() - 4, 4) == ".cif") {
            Cif cif(file_name);
            (*this) = BasicPdb<ModelType>(cif);
        } else {
            die("Please give me a file ended with '.pdb' or '.cif'!");
        }
    }

    ModelType &operator [](int n) {
        return models[n];
    }

    const ModelType &operator [](int n) const {
        return models[n];
    }

    void print_models() {
        int index = 1;
        std::ofstream out;
        for (auto && model : models) {
            std::string pdb_name = name + "-model-" + boost::lexical_cast<std::string>(index);
            write_pdb(model, pdb_name + ".pdb");
            index++;
        }
    }

    void print(std::ostream &out) {
        int index = 1;
        for (auto && model : models) {
            out << "MODEL " << index << "\n";
            out << model;    
            out << "ENDMDL\n";
            index++;
        }
    }

};

typedef BasicPdb<Model> Pdb;

} /// namespace jian

#endif // BasicPdb<ModelType>_H

