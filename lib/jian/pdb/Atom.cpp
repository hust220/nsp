#include "Atom.hpp"

namespace jian {

std::map<char, double> Atom::s_mass {
    {'H', 1.00794},
    {'C', 12.0107},
    {'O', 15.9994},
    {'N', 14.0067},
    {'P', 30.973762},
    {'S', 32},
};
 
Atom::Atom(std::string name, double x, double y, double z) {
    set_name(name);
    at(0) = x; at(1) = y; at(2) = z;
}

Atom::Atom(const std::string &name, int num, double x, double y, double z) {
    set_name(name);
    this->num = num;
    at(0) = x; at(1) = y; at(2) = z;
}

void Atom::set_name(const std::string &s) {
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
//        std::cout << "Unknown atom: " << name[0] << std::endl;
//        throw e;
    }
}

} // namespace jian

