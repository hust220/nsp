#pragma once

#include "pdb_molecule.hpp"
#include "pdb_parser.hpp"
#include "file.hpp"

BEGIN_JN

class MolWriter {
public:
    const static std::vector<std::string> chain_names;
    int atom_num, residue_num, model_num;
    S atom_name, residue_name, chain_name, atom_name_label;
    double x, y, z, a, b;
    std::ostream stream;

    MolWriter();

    MolWriter(std::ostream &out);

    void init();

    void bind_stream(std::ostream &s);

    void write();

    void read(const Atom &atom);

    void write(const Atom &atom);

    void write(const Residue &residue);

    void write(const Chain &chain);

    template<typename F>
    void write_model(F &&f) {
        write_model_begin();
        f();
        write_model_end();
        model_num++;
        atom_num = 1;
        residue_name = "X";
    }

    void write(const Model &model);

    void write(const Molecule &mol);

    void write_model_begin();

    void write_model_end();

    void write_file_end();

    void write_chain_end();

};

void chain_read_model(Chain &chain, S f, S type = "");
Chain read_model_to_chain(S f, S type = "");
//void append_chain_to_file(const Chain &chain, const S &file_name, int n);

template<typename T>
void mol_write(const T & t, S file_name) {
    std::ofstream(file_name.c_str()) << t;
}

template<typename T>
void mol_read(T &t, S file_name, S mol_type = "") {
    MolParser *parser = MolParser::make(file::type(file_name), file_name, mol_type);
    (*parser) >> t;
    delete parser;
}

template<typename T>
T mol_read_to(S f, S type = "") {
    T t;
    mol_read(t, f, type);
    return t;
}

MolParser &operator >> (MolParser &input, Atom &atom);
MolParser &operator >> (MolParser &input, Residue &residue);
MolParser &operator >> (MolParser &input, Chain &chain);
MolParser &operator >> (MolParser &input, Model &model);
MolParser &operator >> (MolParser &input, Molecule &mol);

std::ostream &operator <<(std::ostream &output, const Atom &atom);
std::ostream &operator <<(std::ostream &output, const Residue &residue);
std::ostream &operator <<(std::ostream &output, const Chain &chain);
std::ostream &operator <<(std::ostream &output, const Model &model);
std::ostream &operator <<(std::ostream &output, const Molecule &mol);

class MolReader {
public:
    class ModelEl : public BasicEl<ModelEl, Model> {
    public:
        using Data = Model;
        using El = ModelEl;
    };

    class ModelIt : public BasicIt<ModelIt, ModelEl>
    {
    public:
        using Data = Model;
        using El = ModelEl;
        using It = ModelIt;

        SP<El> el;
        int n;
        MolParser *parser;

        ModelIt() : 
            n(-1), parser(NULL)
        {}

        ModelIt(MolParser *parser_) :
            n(-1), parser(parser_)
        {
            (*this)++;
        }

        virtual Data &operator *() const
        {
            return el->data;
        }

        virtual bool operator ==(It other) const
        {
            return n == other.n;
        }

        It &operator ++()
        {
            if (parser == NULL || parser->eof()) throw "";
            el = STD_ make_shared<El>();
            (*parser) >> el->data;
            n = (parser->eof() ? -1 : n + 1);
            return *this;
        }

    };

    SP<MolParser> parser;

    MolReader(Str filename, Str mol_type = "") {
        parser.reset(MolParser::make(jian::file::type(filename), filename, mol_type));
    }

    ModelIt model_begin() const {
        return ModelIt(parser.get());
    }

    ModelIt model_end() const {
        return ModelIt{};
    }
};

template<typename F>
void for_each_model(S filename, F && f, S mol_type = "") {
    int i = 0;
    MolParser *parser = MolParser::make(jian::file::type(filename), filename, mol_type);
    Model m;

    do {
        (*parser) >> m;
        f(m, i);
        i++;
        m.clear();
    } while (!parser->eof());

    delete parser;
}

template<typename T, typename U, typename = std::enable_if_t<!(std::is_same<T, Residue>::value), int>>
T coarse_grained(const T & t, const U & ls) {
    T t_;
    t_.name = t.name;
    for (auto && r : t) {
        t_.push_back(coarse_grained(r, ls));
    }
    return t_;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_same<T, Residue>::value, int>>
T coarse_grained(const T &residue, U && ls) {
    Residue r;
    r.name = residue.name;
    r.num = residue.num;
    for (auto && atom : residue) {
        if (std::find(ls.begin(), ls.end(), atom.name) != ls.end()) {
            r.push_back(atom);
        }
    }
    return r;
}

END_JN

