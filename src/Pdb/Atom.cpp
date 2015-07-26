#include "Atom.h"

namespace jian {

Atom::Atom(string &line, string rnaName) {
    // set name
    name = line.substr(12, 4);
    string str;
    for (int i = 0; i < (int) name.size(); i++) {
        if (name[i] == ' ') {
            continue;
        } else if (name[i] == '\'') {
            str += '*';
        } else {
            str += name[i];
        }
    }
    name = str;
    if (name == "OP1") {
        name = "O1P";
    } else if (name == "OP2") {
        name = "O2P";
    }

    /* set resName */
    resName = line.substr(17, 3);
    str = "";
    for (int i = 0; i < (int) resName.size(); i++) {
        if (resName[i] != ' ') {
            str += resName[i];
        }
    }
    resName = str;

    /* set num */
    num = atoi(line.substr(6, 5).c_str());

    /* set rnaName */
    this->rnaName = rnaName;

    /* set line */
    this->line = line;

    /* set mass */
    if (name[0] == 'H') {
        mass = MASS_H;
    }    else if (name[0] == 'C') {
        mass = MASS_C;
    }    else if (name[0] == 'O') {
        mass = MASS_O;
    }    else if (name[0] == 'N') {
        mass = MASS_N;
    }    else if (name[0] == 'P') {
        mass = MASS_P;
    }

    // set atom's coordinates
    x = atof(line.substr(30, 8).c_str());
    y = atof(line.substr(38, 8).c_str());
    z = atof(line.substr(46, 8).c_str());
}

Atom::Atom(const string &name, double x, double y, double z) {
    this->name = name;
    replace(this->name.begin(), this->name.end(), '\'', '*');
    if (name == "OP1") {
        this->name = "O1P";
    } else if (name == "OP2") {
        this->name = "O2P";
    }

    this->x = x;
    this->y = y;
    this->z = z;

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

Atom::Atom(Point p, string resName, string name, int num) {
    this->name = name;
    replace(this->name.begin(), this->name.end(), '\'', '*');
    if (name == "OP1") {
        this->name = "O1P";
    } else if (name == "OP2") {
        this->name = "O2P";
    }

    this->resName = resName;
    this->num = num;
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
    x = p.x;
    y = p.y;
    z = p.z;
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


