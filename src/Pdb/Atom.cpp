#include "Atom.h"

namespace jian {

Atom::Atom() {
}

Atom::Atom(MolFile &pdb_file) {
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

Atom::Atom(string &line) {
    set_name(line.substr(12, 4));
    num = atoi(line.substr(6, 5).c_str());
    set_mass();
    x = atof(line.substr(30, 8).c_str());
    y = atof(line.substr(38, 8).c_str());
    z = atof(line.substr(46, 8).c_str());
}

Atom::Atom(std::string name, double x, double y, double z) {
    set_name(name);
    this->x = x;
    this->y = y;
    this->z = z;
    set_mass();

}

Atom::Atom(Point p, string resName, string name, int num) {
    set_name(name);
    this->num = num;
    set_mass();
    x = p.x;
    y = p.y;
    z = p.z;
}

void Atom::set_name(std::string name) {
    this->name = name;
    replace(this->name.begin(), this->name.end(), '\'', '*');
    if (name == "OP1") {
        this->name = "O1P";
    } else if (name == "OP2") {
        this->name = "O2P";
    }
}

void Atom::set_mass() {
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

ostream &operator <<(ostream &output, const Atom &atom) {
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

Point *Atom::coord() {
    Point *p = new Point(x, y, z);
    return p;
}

double Atom::dist(Atom &a) {
    double dx = x - a.x;
    double dy = y - a.y;
    double dz = z - a.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

} /// namespace jian


