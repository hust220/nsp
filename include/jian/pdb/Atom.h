#ifndef JIAN_PDB_ATOM_H
#define JIAN_PDB_ATOM_H

#include "../etl.h"
#include "PdbFile.h"
#include "Cif.h"

namespace jian {

class Atom {
public:
    string name;
    string resName;
    string line;
    int num;
    double mass;
    double x, y, z;
    int flag;

    Atom() {}

    Atom(MolFile &pdb_file) {
        if (!pdb_file.eof()) {
            set_name(pdb_file.atom_name());
            num = pdb_file.atom_num();
            set_mass();
            x = pdb_file.x();
            y = pdb_file.y();
            z = pdb_file.z();
            pdb_file.next();
        }
    }

    Atom(Atom *atom) : name(atom->name), resName(atom->resName),
                            line(atom->line), num(atom->num), mass(atom->mass), 
                            x(atom->x), y(atom->y), z(atom->z), flag(atom->flag) {}

    Atom(const Atom &atom) : name(atom.name), resName(atom.resName), 
                            line(atom.line), num(atom.num), mass(atom.mass), 
                            x(atom.x), y(atom.y), z(atom.z), flag(atom.flag) {}

    Atom &operator =(const Atom &atom) {
        name = atom.name;
        resName = atom.resName;
        line = atom.line;
        num = atom.num;
        mass = atom.mass;
        x = atom.x;
        y = atom.y;
        z = atom.z;
        flag = atom.flag;
    }

    Atom(const Point &point, const string &name) {
        x = point.x;
        y = point.y;
        z = point.z;
        this->name = name;
    }

    friend ostream &operator <<(ostream &, const Atom &);

    double &operator [](int i) {
        if (i == 0) {
            return x;
        } else if (i == 1) {
            return y;
        } else if (i == 2) {
            return z;
        }
    }

    const double &operator [](int i) const {
        if (i == 0) {
            return x;
        } else if (i == 1) {
            return y;
        } else if (i == 2) {
            return z;
        }
    }

    template<class T = Point> T pos() const {
        T t;
        t[0] = x;
        t[1] = y;
        t[2] = z;
        return t;
    }

    Atom(const string &line) {
        set_name(line.substr(12, 4));
        num = atoi(line.substr(6, 5).c_str());
        set_mass();
        x = atof(line.substr(30, 8).c_str());
        y = atof(line.substr(38, 8).c_str());
        z = atof(line.substr(46, 8).c_str());
    }

    Atom(std::string name, double x, double y, double z) {
        set_name(name);
        this->x = x;
        this->y = y;
        this->z = z;
        set_mass();

    }

    Atom(Point p, string resName, string name, int num) {
        set_name(name);
        this->num = num;
        set_mass();
        x = p.x;
        y = p.y;
        z = p.z;
    }

    void set_name(std::string name) {
        this->name = name;
        replace(this->name.begin(), this->name.end(), '\'', '*');
        if (name == "OP1") {
            this->name = "O1P";
        } else if (name == "OP2") {
            this->name = "O2P";
        }
    }

    void set_mass() {
        if (name[0] == 'H') {
            mass = MASS_H;
        } else if (name[0] == 'C') {
            mass = MASS_C;
        } else if (name[0] == 'O') {
            mass = MASS_O;
        } else if (name[0] == 'N') {
            mass = MASS_N;
        } else if (name[0] == 'P') {
            mass = MASS_P;
        } else {
        }
    }

    Point *coord() {
        Point *p = new Point(x, y, z);
        return p;
    }

    double dist(Atom &a) {
        double dx = x - a.x;
        double dy = y - a.y;
        double dz = z - a.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    int empty() {
        return name == "";
    }

};

inline ostream &operator <<(ostream &output, const Atom &atom) {
    output << fixed << setprecision(3);
    output << "ATOM" 
           << setw(7)  << 1 << "  "
           << left << setw(4)  << atom.name
           << right << setw(3) << "A"
           << setw(2)  << "A"
           << setw(4)  << 1 
           << setw(12) << atom.x 
           << setw(8)  << atom.y 
           << setw(8)  << atom.z 
           << "\n";
    return output;
}

inline Atom make_atom(const std::string &line) {
    Atom atom;
    atom.set_name(line.substr(12, 4));
    atom.num = boost::lexical_cast<int>(line.substr(6, 5));
    atom.set_mass();
    atom.x = boost::lexical_cast<double>(line.substr(30, 8));
    atom.y = boost::lexical_cast<double>(line.substr(38, 8));
    atom.z = boost::lexical_cast<double>(line.substr(46, 8));
    return atom;
}

} // namespace jian

#endif
