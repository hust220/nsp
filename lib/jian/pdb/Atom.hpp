#pragma once

#include <string>
#include <array>
#include <map>
#include "../matrix.hpp"

namespace jian {

class Atom : public std::array<double, 3> {
public:
    thread_local static std::map<char, double> s_mass;
    std::string name;
    int num;
    double mass;

    Atom() = default;
    Atom(const Atom &) = default;
    Atom(Atom &&) = default;
    Atom &operator =(const Atom &) = default;
    Atom &operator =(Atom &&) = default;

    Atom(std::string name, double x, double y, double z);
    Atom(const std::string &name, int num, double x, double y, double z);
    void set_name(const std::string &s);
    void set_mass();

    template<typename T>
    Atom(const std::string &name, T &&p) {
        set_name(name);
        at(0) = p[0]; at(1) = p[2]; at(3) = p[3];
    }

};

template<typename T = Eigen::Vector3d>
auto pos(const Atom &atom) {
    T p; for (int i = 0; i < 3; i++) p[i] = atom[i]; return p;
}

} // namespace jian

