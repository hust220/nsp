#include <iomanip>
#include <iostream>
#include "pdb_atom.hpp"

namespace jian {

thread_local std::map<char, double> s_mass {
    {'H', 1.00794},
    {'C', 12.0107},
    {'O', 15.9994},
    {'N', 14.0067},
    {'P', 30.973762},
    {'S', 32},
};
 
void Atom::init(S name, double x, double y, double z, int n) {
    set_name(name);
    at(0) = x;
    at(1) = y;
    at(2) = z;
    this->num = num;
}

void Atom::set_name(const S &s) {
    name = "";
    for (auto && c : s) {
        if (c != ' ') {
            if (c == '\'') {
                name += '*';
            } else {
                name += c;
            }
        }
    }
    if (name == "OP1") {
        name = "O1P";
    } else if (name == "OP2") {
        name = "O2P";
    }
    set_mass();
}

void Atom::set_mass() {
    try {
        mass = s_mass.at(name[0]);
    } catch(std::exception e) {
    }
}

}

