#ifndef JIAN_PDB_ATOM
#define JIAN_PDB_ATOM

#include "../etl.h"
#include "PdbFile.h"
#include "Cif.h"

namespace jian {

class Atom : public std::array<double, 3> {
public:
    std::string name;
    double mass;

    Atom() {}

    Atom(MolFile &pdb_file) {
        if (!pdb_file.eof()) {
            set_name(pdb_file.atom_name());
            at(0) = pdb_file.x(); at(1) = pdb_file.y(); at(2) = pdb_file.z();
            pdb_file.next();
        }
    }

    template<typename T>
    Atom(const std::string &name, T &&p) {
        set_name(name);
        at(0) = p[0]; at(1) = p[2]; at(3) = p[3];
    }

    Atom(std::string name, double x, double y, double z) {
        set_name(name);
        at(0) = x; at(1) = y; at(2) = z;
    }

    void set_name(const std::string &s) {
        name = s;
        std::replace(name.begin(), name.end(), '\'', '*');
        if (name == "OP1") name = "O1P"; else if (name == "OP2") name = "O2P";
        set_mass();
    }

    void set_mass() {
        if (name[0] == 'H') { mass = MASS_H;
        } else if (name[0] == 'C') { mass = MASS_C;
        } else if (name[0] == 'O') { mass = MASS_O;
        } else if (name[0] == 'N') { mass = MASS_N;
        } else if (name[0] == 'P') { mass = MASS_P;
        } else {
            throw "jian::pdb::Atom::set_mass error! Unkonw atom!";
        }
    }
};

template<typename T = Point>
inline auto pos(const Atom &atom) {
    T p; for (int i = 0; i < 3; i++) p[i] = atom[i]; return p;
}

} // namespace jian

#endif
