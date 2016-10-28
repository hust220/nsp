#pragma once

#include <string>
#include <array>
#include <map>
#include "../matrix.hpp"

namespace jian {

class Atom : public std::array<double, 3> {
public:
    std::string name;
    int num;
    double mass;

    void init(std::string name, double x, double y, double z, int num = -1);

    void set_name(const std::string &s);
    void set_mass();

};

} // namespace jian

